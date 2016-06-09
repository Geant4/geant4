#!/usr/bin/python3
# ==================================================================
# An example for reading a GDML file
#
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
    self.world= None
    self.gdml_parser= G4GDMLParser()

  # -----------------------------------------------------------------
  def __del__(self):
    pass
    
  # -----------------------------------------------------------------
  def Construct(self):
    self.gdml_parser.Read("qgeom.gdml")
    self.world= self.gdml_parser.GetWorldVolume()

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
gApplyUICommand("/vis/open OGLIX")
gApplyUICommand("/vis/scene/create")
gApplyUICommand("/vis/scene/add/volume")
gApplyUICommand("/vis/sceneHandler/attach")
gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 90. -90.")

