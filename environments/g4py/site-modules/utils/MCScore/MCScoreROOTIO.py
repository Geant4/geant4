"""
Python module

This module provides ROOT IO interface for MCScore data

  [C] MCScoreROOTIO
  [f] loop_tree(tfile, analyze_vertex):

                                              Q, 2006
"""
from array import array
import string
from Geant4.hepunit import *
import mcscore
import ROOT

# ==================================================================
# public symbols
# ==================================================================
__all__ = [ 'MCScoreROOTIO', 'loop_tree' ]


# ==================================================================
#   class definition
# ==================================================================
# ------------------------------------------------------------------
# MCScoreROOTIO
# ------------------------------------------------------------------
class MCScoreROOTIO:
  "ROOT IO interface for MCScore"
  def __init__(self, bsize=250):
    self.maxparticle = bsize # buffer size for #particles/vertex

  # ------------------------------------------------------------------
  def define_tree(self):
    "define ROOT tree"
    # defining tree header...
    self.vtree = ROOT.TTree('vertex',   'mc vertex')
    self.ptree = ROOT.TTree('particle', 'mc particle')

    # vertex...  
    self.a_x = array('d', [0.]); self.vtree.Branch('x', self.a_x, 'x/d')
    self.a_y = array('d', [0.]); self.vtree.Branch('y', self.a_y, 'y/d')
    self.a_z = array('d', [0.]); self.vtree.Branch('z', self.a_z, 'z/d')
    
    #global a_np, a_namelist, a_Z, a_A, a_ke, a_px, a_py, a_pz
    self.a_np = array('i', [0]); self.ptree.Branch('np', self.a_np, 'np/i')
    
    self.a_namelist = array('c', self.maxparticle*10*['\0'])
      # 10 characters/particle
    self.ptree.Branch('namelist', self.a_namelist, 'namelist/C')
    
    self.a_Z = array('i', self.maxparticle*[0])
    self.ptree.Branch('Z',  self.a_Z,  'Z[np]/I')

    self.a_A = array('i', self.maxparticle*[0])
    self.ptree.Branch('A',  self.a_A,  'A[np]/i')

    self.a_ke = array('d', self.maxparticle*[0.])
    self.ptree.Branch('kE', self.a_ke, 'kE[np]/d')

    self.a_px = array('d', self.maxparticle*[0.])
    self.ptree.Branch('px', self.a_px, 'px[np]/d')

    self.a_py = array('d', self.maxparticle*[0.])
    self.ptree.Branch('py', self.a_py, 'py[np]/d')

    self.a_pz = array('d', self.maxparticle*[0.])
    self.ptree.Branch('pz', self.a_pz, 'pz[np]/d')

  # ------------------------------------------------------------------
  def fill_tree(self, vertex):
    "fill vertex information to ROOT tree"
    # ------------------------------------------------------------------
    def push_pname(i0, pname): # local function
      n = len(pname)
      for i in xrange(n):
        self.a_namelist[i0+i] = pname[i]
        self.a_namelist[i0+n] = ' '

    self.a_x[0] = vertex.x
    self.a_y[0] = vertex.y
    self.a_z[0] = vertex.z
    self.vtree.Fill()
          
    if vertex.nparticle > self.maxparticle:
      raise """
      *** buffer overflow in #particles/vertex.
      *** please increment buffersize in MCScoreROOTIO(bsize).
      """, self.maxparticle
    
    idx_namelist = 0
    self.a_np[0] = vertex.nparticle
    for ip in range(vertex.nparticle):
      particle = vertex.particle_list[ip]
      push_pname(idx_namelist, particle.name)
      idx_namelist += (len(particle.name)+1)
      self.a_Z[ip] = particle.Z
      self.a_A[ip] = particle.A
      self.a_ke[ip] = particle.kineticE
      self.a_px[ip] = particle.px
      self.a_py[ip] = particle.py
      self.a_pz[ip] = particle.pz

    self.a_namelist[idx_namelist] = '\0'
    self.ptree.Fill()
      
                  
# ==================================================================
# functions
# ==================================================================
def loop_tree(tfile, analyze_vertex):
  """
  loop ROOT tree in a ROOT file.
    * analyze_vertex: user function : analyze_vertex(MCVertex)
  """
  avtree = tfile.Get("vertex")
  aptree = tfile.Get("particle")

  # reading vertex...
  n_vertex = avtree.GetEntries()
  for ivtx in xrange(n_vertex):
    avtree.GetEntry(ivtx)
    aptree.GetEntry(ivtx)

    # vertex
    vertex = MCScore.MCVertex(avtree.x, avtree.y, avtree.z)    

    # reading secondary particles...
    nsec = aptree.np
    namelist = aptree.namelist
    pname = namelist.split()
    for ip in xrange(nsec):
      particle = MCScore.MCParticle(pname[ip], aptree.Z[ip], aptree.A[ip],
                                    aptree.kE[ip], aptree.px[ip],
                                    aptree.py[ip], aptree.pz[ip])
      vertex.append_particle(particle)

    analyze_vertex(vertex)

  return n_vertex


# ==================================================================
# test
# ==================================================================
def test_t2root():
  g = ROOT.TFile("reaction.root", 'recreate')

  rootio = MCScoreROOTIO()
  rootio.define_tree()

  f = open("reaction.dat")
  f.seek(0)

  while(1):
    vertex = MCScore.read_next_vertex(f)
    if vertex == 0:
      break

    # filling ...
    rootio.fill_tree(vertex)

    del vertex

  print ">>> EOF"
  f.close()  

  # closing ROOT file
  g.Write()
  g.Close()

# ------------------------------------------------------------------
def test_loop():
  def my_analysis(vertex):
    vertex.printout()

  f = ROOT.TFile("reaction.root", 'read')  
  nv = loop_tree(f, my_analysis)
  print "*** # of vertex= ", nv
  

# ==================================================================
# main
# ==================================================================
if __name__ == "__main__":
  test_t2root()
  test_loop()

