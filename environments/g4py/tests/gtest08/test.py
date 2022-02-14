#!/usr/bin/python3
# ==================================================================
# python script for Geant4Py test
#
#   gtest09
#   - test for checking use of G4PhysListFactory
# ==================================================================

from Geant4 import *

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

    return pv_world

class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
  "My Primary Generator Action"

  def __init__(self):
    G4VUserPrimaryGeneratorAction.__init__(self)
    self.particleGun= G4ParticleGun(1)

  def GeneratePrimaries(self, event):
    dx=0.
    self.particleGun.SetParticleMomentumDirection(G4ThreeVector(dx, 0., 1.))
    self.particleGun.GeneratePrimaryVertex(event)

# ==================================================================
# main
# ==================================================================
# set geometry
myDC = MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)

# Create physics list
factory = G4PhysListFactory()
myPhysList = factory.GetReferencePhysList("FTFP_BERT")

if myPhysList is None:
  raise RuntimeError("No physics list named FTFP_BERT found")

gRunManager.SetUserInitialization(myPhysList)

# Event generator
myGenAction = MyPrimaryGeneratorAction()
myGenAction.particleGun.SetParticleByName("e-")
myGenAction.particleGun.SetParticleEnergy(200.*MeV)
myGenAction.particleGun.SetParticlePosition(G4ThreeVector(0.,0.,-14.9)*cm)

gRunManager.SetUserAction(myGenAction)

# Init, run terminate
gRunManager.Initialize()
gRunManager.BeamOn(10)
gTerminate()
