#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
# ==================================================================
from Geant4 import *
import demo_wp

# ==================================================================
# user actions in python
# ==================================================================
class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
  "My Primary Generator Action"

  def __init__(self):
    G4VUserPrimaryGeneratorAction.__init__(self)
    self.particleGun= G4ParticleGun(1)

  def GeneratePrimaries(self, event):
    self.particleGun.GeneratePrimaryVertex(event)

# ------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
  "My Run Action"
  
  def EndOfRunAction(self, run):
    print "*** End of Run"
    print "- Run sammary : (id= %d, #events= %d)" \
          % (run.GetRunID(), run.GetNumberOfEventToBeProcessed())

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
  "My Event Action"

  def EndOfEventAction(self, event):
    pass

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
  "My Stepping Action"

  def UserSteppingAction(self, step):
    pass
    #print "*** dE/dx in current step=", step.GetTotalEnergyDeposit()
    preStepPoint= step.GetPreStepPoint()
    track= step.GetTrack()
    touchable= track.GetTouchable()
    #print "*** vid= ", touchable.GetReplicaNumber()


# ==================================================================
# main
# ==================================================================
myMaterials= demo_wp.MyMaterials()
myMaterials.Construct()

myDC= demo_wp.MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)

myPL= demo_wp.MyPhysicsList()
gRunManager.SetUserInitialization(myPL)

# set user actions...
myPGA= MyPrimaryGeneratorAction()
gRunManager.SetUserAction(myPGA)

myRA= MyRunAction()
gRunManager.SetUserAction(myRA)
  
myEA= MyEventAction()
gRunManager.SetUserAction(myEA)

mySA= MySteppingAction()
gRunManager.SetUserAction(mySA)


# set particle gun
pg= myPGA.particleGun
pg.SetParticleByName("proton")
pg.SetParticleEnergy(230.*MeV)
pg.SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.))
pg.SetParticlePosition(G4ThreeVector(0.,0.,-20.)*cm)
  
gRunManager.Initialize()

# visualization
gApplyUICommand("/control/execute vis.mac")
  
# beamOn
#gRunManager.BeamOn(3)


