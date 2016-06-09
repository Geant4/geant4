"""
Python module

This module provides classes and functions for scoring reactions

  [C] MCVertex:
  [C] MCParticle:
  [f] read_next_vertex(stream):

                                              Q, 2006
"""
import string
from Geant4.hepunit import *

# ==================================================================
# public symbols
# ==================================================================
__all__ = [ 'MCParticle', 'MCVertex', 'read_next_vertex' ]


# ==================================================================
#   class definition
# ==================================================================

# ------------------------------------------------------------------
# MCParticle
# ------------------------------------------------------------------
class MCParticle:
  "MC particle"
  def __init__(self, aname, aZ, aA, akE, apx, apy, apz):
    self.name = aname
    self.Z = aZ
    self.A = aA
    self.kineticE = akE
    self.px = apx
    self.py = apy
    self.pz = apz

  def printout(self):
    print "--- particle: %s, Z=%2d, A=%2d, kE=%g" % \
          (self.name, self.Z, self.A, self.kineticE/MeV)

# ------------------------------------------------------------------
# MCVertex
# ------------------------------------------------------------------
class MCVertex :
  "MC vertex"
  def __init__(self, ax, ay, az):
    self.x = ax
    self.y = ay
    self.z = az
    self.nparticle = 0
    self.particle_list = []

  def append_particle(self, aparticle):
    self.particle_list.append(aparticle)
    self.nparticle= self.nparticle+1

  def printout(self):
    print "@@@ vertex: x=(%g,%g,%g) Nsec=%3d" % \
          (self.x/cm, self.y/cm, self.z/cm, self.nparticle)
    for p in self.particle_list:
      p.printout()
  
  def dump_vertex(self, stream):
    aline = "%g %g %g %d\n" % \
            (self.x/m, self.y/m, self.z/m, self.nparticle)
    stream.write(aline)
    for p in self.particle_list:
      aline = " %s %d %d %g %g %g %g\n" % \
              (p.name, p.Z, p.A, p.kineticE/MeV, p.px/MeV, p.py/MeV, p.pz/MeV)
      stream.write(aline)

  def __del__(self):
    np = len(self.particle_list)
    del self.particle_list[0:np]
    

# ==================================================================
#   I/O interface
# ==================================================================
def read_next_vertex(stream):
  "read next vertex from a file stream"
  line= stream.readline()
  if line == "":  # EOF
    return 0

  # reading vertex
  data = line.split()
  x = string.atof(data[0]) * m
  y = string.atof(data[1]) * m
  z = string.atof(data[2]) * m
  nsec = string.atoi(data[3])

  vertex = MCVertex(x,y,z)

  # reading particles
  for p in range(0, nsec):
    data = stream.readline().split()
    pname = data[0]
    Z = string.atoi(data[1])
    A = string.atoi(data[2])    
    kE = string.atof(data[3]) * MeV
    px = string.atof(data[4]) * MeV
    py = string.atof(data[5]) * MeV
    pz = string.atof(data[6]) * MeV

    particle = MCParticle(pname, Z, A, kE, px, py, pz)
    vertex.append_particle(particle)

  return vertex


# ==================================================================
# test
# ==================================================================
def test():
  f = open("reaction.dat")
  f.seek(0)

  while(1):
    vertex = read_next_vertex(f)
    if vertex == 0:
      break
    vertex.printout()
    del vertex
  f.close()
  print ">>> EOF"


# ==================================================================
# main
# ==================================================================
if __name__ == "__main__":
  test()
  
