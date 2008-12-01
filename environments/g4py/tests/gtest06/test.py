#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest06
#   - test for constructing/visualizing boolean geoemtries
# ==================================================================
from Geant4 import *
import g4py.ExN01pl, g4py.ParticleGun

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

    self.sld_box= G4Box("box",20.*cm, 20.*cm, 20.*cm);
    self.sld_cyl= G4Tubs("cylinder",0., 10.*cm, 30.*cm, 0., twopi)

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
  def ConstructUnion(self):
    global sld_union
    sld_union= G4UnionSolid("box+cylinder", self.sld_box, self.sld_cyl); 
 
    self.lv_object.SetSolid(sld_union)
    self.lv_object.SetVisAttributes(self.va_blue)
    gRunManager.GeometryHasBeenModified()
    
  # -----------------------------------------------------------------
  def ConstructIntersection(self):
    offset= G4ThreeVector(20.*cm, 20.*cm, 0.)
    global sld_intersection  
    sld_intersection= G4IntersectionSolid("box*cylinder",
                                          self.sld_box, self.sld_cyl,
                                          None, offset)

    self.lv_object.SetSolid(sld_intersection)
    self.lv_object.SetVisAttributes(self.va_magenta)
    gRunManager.GeometryHasBeenModified()

  # -----------------------------------------------------------------
  def ConstructSubtraction(self):
    global sld_subtraction  
    sld_subtraction= G4SubtractionSolid("box-cylinder",
                                        self.sld_box, self.sld_cyl)

    self.lv_object.SetSolid(sld_subtraction)
    self.lv_object.SetVisAttributes(self.va_red)
    gRunManager.GeometryHasBeenModified()
            
  # -----------------------------------------------------------------
  def Construct(self): # return the world volume
    return self.world  
  
# ==================================================================
# main
# ==================================================================
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
gApplyUICommand("/vis/open RayTracer")
gApplyUICommand("/vis/rayTracer/headAngle 40.")
gApplyUICommand("/vis/rayTracer/eyePosition 100 100 150  cm")

# create a vrml file for each solid type
f_list= (
  ("union",          myDC.ConstructUnion),
  ("intersection",   myDC.ConstructIntersection),
  ("subtraction",    myDC.ConstructSubtraction)
  )


for s,f in f_list:
  f.__call__()
  fname= "%s.jpg" % (s)
  cmdstr= "/vis/rayTracer/trace " + fname
  gApplyUICommand(cmdstr)
  

