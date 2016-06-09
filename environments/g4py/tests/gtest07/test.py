#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest07
#   - test for checking overlapped geometries
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

  # -----------------------------------------------------------------
  def Construct(self):
    # Python has automatic garbage collection system.
    # Geometry objects must be defined as GLOBAL not to be deleted.
    air= gNistManager.FindOrBuildMaterial("G4_AIR")

    # world volume
    global sld_world, lv_world, pv_world
    sld_world= G4Box("world", 1.*m, 1.*m, 1.*m)
    lv_world= G4LogicalVolume(sld_world, air, "world")    
    pv_world= G4PVPlacement(G4Transform3D(), lv_world, "world",
                            None, False, 0)

    # box
    global sld_box, lv_box, pv_box
    sld_box= G4Box("box", 10.*cm, 10.*cm, 10.*cm);
    lv_box= G4LogicalVolume(sld_box, air, "box")
    pv_box= G4PVPlacement(None, G4ThreeVector(), "box", lv_box,
                          pv_world, False, 0, True)

    # cylinder
    global sld_cyl, lv_cyl, pv_cyl1, pv_cyl2, pv_cyl3
    sld_cyl= G4Tubs("cylinder",0., 2.*cm, 2.*cm, 0., twopi)
    lv_cyl= G4LogicalVolume(sld_cyl, air, "cylinder")

    #
    # the following placements are !! overlapped !!
    #
    # doubly placed
    pv_cyl1= G4PVPlacement(None, G4ThreeVector(), "cylinder", lv_cyl,
                           pv_world, False, 0, True)

    # overlaped
    pv_cyl2= G4PVPlacement(None, G4ThreeVector(10.*cm,0.,0.),
                           "cylinder", lv_cyl,
                           pv_world, False, 1, True)

    # sticked out
    pv_cyl3= G4PVPlacement(None, G4ThreeVector(10.*cm,0.,0.),
                           "cylinder", lv_cyl,
                           pv_box, False, 0, True)

    return pv_world
  
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
gRunManager.Initialize() # overlap should be detected !!


