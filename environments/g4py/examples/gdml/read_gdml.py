#!/usr/bin/python
# ==================================================================
# An example for reading a GDML file
#
# ==================================================================
from Geant4 import *
from Geant4.G4gdml import *

import ExN01pl, ParticleGun

# ==================================================================
# user actions in python
# ==================================================================
class MyDetectorConstruction(G4VUserDetectorConstruction):
  "My Detector Construction"

  def __init__(self):
    G4VUserDetectorConstruction.__init__(self)
    self.world= None
    self.sxp= SAXProcessor()
    self.config= ProcessingConfigurator()

  # -----------------------------------------------------------------
  def __del__(self):
    self.sxp.Finalize()
    
  # -----------------------------------------------------------------
  def Construct(self):
    self.sxp.Initialize()

    self.config.SetURI("qgeom.gdml")
    self.config.SetSetupName("Default")

    self.sxp.Configure(self.config)
    self.sxp.Run()

    self.world= GDMLProcessor.GetInstance().GetWorldVolume()

    return self.world


# ==================================================================
# main
# ==================================================================
# set geometry
myDC= MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)

# minimal physics list
ExN01pl.Construct()

# set primary generator action
ParticleGun.Construct()

# initialize
gRunManager.Initialize()

# visualization
gApplyUICommand("/vis/open OGLIX")
gApplyUICommand("/vis/scene/create")
gApplyUICommand("/vis/scene/add/volume")
gApplyUICommand("/vis/sceneHandler/attach")
gApplyUICommand("/vis/viewer/set/viewpointThetaPhi 90. -90.")

