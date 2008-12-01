#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest05
#   - test for CSG geometry construction in Python
# ==================================================================
from Geant4 import *
import g4py.ExN01pl, g4py.ParticleGun
import os, math

# ==================================================================
# user actions in python
# ==================================================================
class MyDetectorConstruction(G4VUserDetectorConstruction):
  "My Detector Construction"

  def __init__(self):
    G4VUserDetectorConstruction.__init__(self)
    self.air= gNistManager.FindOrBuildMaterial("G4_AIR")
    self.lv_object= None
    self.world= self.ConstructWorld()
    
    self.va_red= G4VisAttributes(G4Color(1.,0.,0.))
    self.va_cyan= G4VisAttributes(G4Color(0.,1.,1.))
    self.va_green= G4VisAttributes(G4Color(0.,1.,0.))
    self.va_blue= G4VisAttributes(G4Color(0.,0.,1.))
    self.va_magenta= G4VisAttributes(G4Color(1.,0.,1.))

  # -----------------------------------------------------------------
  def ConstructWorld(self):
    # Python has automatic garbage collection system.
    # Geometry objects must be defined as GLOBAL not to be deleted.
    global sld_world, lv_world, pv_world, va_world

    sld_world= G4Box("world", 1.*m, 1.*m, 1.*m)
    lv_world= G4LogicalVolume(sld_world, self.air, "world")    
    pv_world= G4PVPlacement(G4Transform3D(), lv_world, "world",
                            None, False, 0)

    va_world= G4VisAttributes()
    va_world.SetVisibility(False)
    lv_world.SetVisAttributes(va_world)

    # solid object (dummy)
    global sld_sld, lv_sld, pv_sld
    sld_sld= G4Box("dummy", 10.*cm, 10.*cm, 10.*cm)
    self.lv_object= lv_sld= G4LogicalVolume(sld_sld, self.air, "dummy")
    pv_sld= G4PVPlacement(None, G4ThreeVector(), "dummy", lv_sld,
                          pv_world, False, 0)

    return pv_world

  # -----------------------------------------------------------------
  def ConstructBox(self):
    global sld_box
    sld_box= G4Box("box", 30.*cm, 40.*cm, 60.*cm)
    self.lv_object.SetSolid(sld_box)
    self.lv_object.SetVisAttributes(self.va_red)
    gRunManager.GeometryHasBeenModified()
    
  # -----------------------------------------------------------------
  def ConstructTubs(self):
    global sld_tubs
    sld_tubs= G4Tubs("tubs", 10.*cm, 15.*cm, 20.*cm, 0., pi)
    self.lv_object.SetSolid(sld_tubs)
    self.lv_object.SetVisAttributes(self.va_cyan)
    gRunManager.GeometryHasBeenModified()
        
  # -----------------------------------------------------------------
  def ConstructCons(self):
    global sld_cons
    sld_cons= G4Cons("cons", 5.*cm, 10.*cm, 20.*cm, 25.*cm,
                     40.*cm, 0., 4./3.*pi)
    self.lv_object.SetSolid(sld_cons)
    self.lv_object.SetVisAttributes(self.va_green)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructPara(self):
    global sld_para
    sld_para= G4Para("para", 30.*cm, 40.*cm, 60.*cm, pi/4., pi/8., 0.)
    self.lv_object.SetSolid(sld_para)
    self.lv_object.SetVisAttributes(self.va_blue)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTrd(self):
    global sld_trd
    sld_trd= G4Trd("trd", 30.*cm, 10.*cm, 40.*cm, 15.*cm, 60.*cm)
    self.lv_object.SetSolid(sld_trd)
    self.lv_object.SetVisAttributes(self.va_blue)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTrap(self):
    global sld_trap
    sld_trap= G4Trap("trap", 60.*cm, 20.*degree, 5.*degree,
                     40.*cm, 30.*cm, 40.*cm, 10.*degree,
                     16.*cm, 10*cm, 14.*cm, 10.*deg)
    self.lv_object.SetSolid(sld_trap)
    self.lv_object.SetVisAttributes(self.va_green)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructSphere(self):
    global sld_sphere
    sld_sphere= G4Sphere("sphere", 100.*cm, 120.*cm, 0., 180.*deg,
                         0., 180.*deg)
    self.lv_object.SetSolid(sld_sphere)
    self.lv_object.SetVisAttributes(self.va_cyan)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructOrb(self):
    global sld_orb
    sld_orb= G4Orb("orb", 100.*cm)    
    self.lv_object.SetSolid(sld_orb)
    self.lv_object.SetVisAttributes(self.va_red)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTorus(self):
    global sld_torus
    sld_torus= G4Torus("torus", 40.*cm, 60.*cm, 200.*cm, 0., 90.*deg)    
    self.lv_object.SetSolid(sld_torus)
    self.lv_object.SetVisAttributes(self.va_magenta)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructPolycone(self):
    zvec= G4doubleVector()
    rinvec= G4doubleVector()
    routvec= G4doubleVector()

    zvec[:]= [ 5.*cm, 7.*cm, 9.*cm, 11.*cm, 25.*cm, 27.*cm, 29.*cm,
               31.*cm, 35.*cm ]
    rinvec[:]= [0.,0.,0.,0.,0.,0.,0.,0.,0.]
    routvec[:]= [ 0., 10.*cm, 10.*cm, 5.*cm, 5.*cm, 10.*cm,
                  10.*cm, 2.*cm, 2.*cm ]

    global sld_pcon
    sld_pcon= CreatePolycone("pcon", 0., twopi, 9, zvec, rinvec,routvec)
    self.lv_object.SetSolid(sld_pcon)
    self.lv_object.SetVisAttributes(self.va_cyan)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructPolyhedra(self):
    zvec= G4doubleVector()
    rinvec= G4doubleVector()
    routvec= G4doubleVector()

    zvec[:]= [ 0., 5.*cm, 8.*cm, 13.*cm, 30.*cm, 32.*cm, 35.*cm ]
    rinvec[:]= [0.,0.,0.,0.,0.,0.,0. ]
    routvec[:]= [ 0., 15.*cm, 15.*cm, 4.*cm, 4.*cm, 10.*cm, 10.*cm ]

    global sld_pgon
    sld_pgon= CreatePolyhedra("pgon", 0., twopi, 5, 7, zvec, rinvec,routvec)
    self.lv_object.SetSolid(sld_pgon)
    self.lv_object.SetVisAttributes(self.va_green)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructEllipticalTube(self):
    global sld_et
    sld_et= G4EllipticalTube("ellipticaltube", 5.*cm, 10.*cm, 20.*cm)    
    self.lv_object.SetSolid(sld_et)
    self.lv_object.SetVisAttributes(self.va_cyan)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructEllipsoid(self):
    global sld_es
    sld_es= G4Ellipsoid("ellipsoid", 10.*cm, 20.*cm, 50.*cm,
                        -10.*cm, 40.*cm)
    self.lv_object.SetSolid(sld_es)
    self.lv_object.SetVisAttributes(self.va_red)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructEllipticalCone(self):
    global sld_ec
    sld_ec= G4EllipticalCone("ellipticalcone", 30.*cm, 60.*cm,
                             50.*cm, 25.*cm)    
    self.lv_object.SetSolid(sld_ec)
    self.lv_object.SetVisAttributes(self.va_magenta)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructHype(self):
    global sld_hype
    sld_hype= G4Hype("hype", 20.*cm, 30.*cm, 0.7, 0.7, 50.*cm)
    self.lv_object.SetSolid(sld_hype)
    self.lv_object.SetVisAttributes(self.va_blue)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTet(self):
    global sld_tet
    p1= G4ThreeVector(0., 0., math.sqrt(3.)*cm)
    p2= G4ThreeVector(0., 2*math.sqrt(2./3.)*cm, -1./math.sqrt(3)*cm)
    p3= G4ThreeVector(-math.sqrt(2.)*cm, -math.sqrt(2./3.)*cm,
                      -1./math.sqrt(3)*cm)
    p4= G4ThreeVector(math.sqrt(2)*cm, -math.sqrt(2./3.)*cm,
                      -1./math.sqrt(3)*cm)

    sld_tet= G4Tet("tet", 20.*p1, 20.*p2, 20.*p3, 20.*p4)
    self.lv_object.SetSolid(sld_tet)
    self.lv_object.SetVisAttributes(self.va_green)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTwistedBox(self):
    global sld_twb
    sld_twb= G4TwistedBox("twistedbox", 30.*deg, 30.*cm, 40.*cm, 60.*cm)
    self.lv_object.SetSolid(sld_twb)
    self.lv_object.SetVisAttributes(self.va_cyan)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTwistedTrap(self):
    global sld_twtrp
    sld_twtrp= G4TwistedTrap("twistedtrap", 30.*deg,
                             60.*cm, 20.*deg, 5.*deg,
                             40.*cm, 30.*cm, 40.*cm,
                             16.*cm, 10.*cm, 14.*cm, 10.*deg)
    self.lv_object.SetSolid(sld_twtrp)
    self.lv_object.SetVisAttributes(self.va_blue)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTwistedTrd(self):
    global sld_twtrd
    sld_twtrd= G4TwistedTrd("twistedtrd", 30.*cm, 10.*cm,
                            40.*cm, 15.*cm, 60.*cm, 30.*deg)
    self.lv_object.SetSolid(sld_twtrd)
    self.lv_object.SetVisAttributes(self.va_green)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructTwistedTubs(self):
    global sld_twt
    sld_twt= G4TwistedTubs("twistedtube", 60.*deg,
                           10.*cm, 15.*cm, 20.*cm, 90.*deg)
    self.lv_object.SetSolid(sld_twt)
    self.lv_object.SetVisAttributes(self.va_magenta)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def Construct(self): # return the world volume
    return self.world  
  
# ==================================================================
# main
# ==================================================================
os.environ["G4VRML_DEST_DIR"]= "."
os.environ["G4VRMLFILE_MAX_FILE_NUM"]= "1"
os.environ["G4VRMLFILE_VIEWER"]= "echo"

# set geometry
myDC= MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)

# minimal physics list
g4py.ExN01pl.Construct()

# set primary generator action
g4py.ParticleGun.Construct()

# initialize
gRunManager.Initialize()

# visualization
gApplyUICommand("/vis/open VRML2FILE")
gApplyUICommand("/vis/scene/create")
gApplyUICommand("/vis/scene/add/volume")
gApplyUICommand("/vis/sceneHandler/attach")
gApplyUICommand("/vis/scene/add/axes 0. 0. 0. 10. cm")

# create a vrml file for each solid type
f_list= (
  ("g4box",            myDC.ConstructBox),
  ("g4tubs",           myDC.ConstructTubs),
  ("g4cons",           myDC.ConstructCons),
  ("g4para",           myDC.ConstructPara),
  ("g4trd",            myDC.ConstructTrd),
  ("g4trap",           myDC.ConstructTrap),
  ("g4sphere",         myDC.ConstructSphere),
  ("g4orb",            myDC.ConstructOrb),
  ("g4torus",          myDC.ConstructTorus),
  ("g4polycone",       myDC.ConstructPolycone),
  ("g4polyhedra",      myDC.ConstructPolyhedra),
  ("g4ellipticaltube", myDC.ConstructEllipticalTube),
  ("g4ellipsoid",      myDC.ConstructEllipsoid),
  ("g4ellipticalcone", myDC.ConstructEllipticalCone),
  ("g4hype",           myDC.ConstructHype),
  ("g4tet",            myDC.ConstructTet),
  ("g4twistedbox",     myDC.ConstructTwistedBox),
  ("g4twistedtrap",    myDC.ConstructTwistedTrap),
  ("g4twistedtrd",     myDC.ConstructTwistedTrd),
  ("g4twistedtubs",    myDC.ConstructTwistedTubs)
  )

for s,f in f_list:
  f.__call__()
  gRunManager.BeamOn(1)  
  fname= "%s.wrl" % (s)
  os.rename("g4_00.wrl", fname)
  
