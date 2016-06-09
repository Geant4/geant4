#!/usr/bin/python3
# ==================================================================
# python script for Geant4Py test
#
#   gtest01
#   - check basic control flow
# ==================================================================
from Geant4 import *
import gtest01
import random

# ==================================================================
# user actions in python
# ==================================================================
class MyPrimaryGeneratorAction(G4VUserPrimaryGeneratorAction):
  "My Primary Generator Action"

  def __init__(self):
    G4VUserPrimaryGeneratorAction.__init__(self)
    self.particleGun= G4ParticleGun(1)

  def GeneratePrimaries(self, event):
    #dx= random.gauss(0., 0.1)
    dx=0.
    self.particleGun.SetParticleMomentumDirection(G4ThreeVector(dx, 0., 1.))
    self.particleGun.GeneratePrimaryVertex(event)

# ------------------------------------------------------------------
class MyRunAction(G4UserRunAction):
  "My Run Action"

  def BeginOfRunAction(self, run):
    print("*** #event to be processed (BRA)=",
    run.GetNumberOfEventToBeProcessed())

  def EndOfRunAction(self, run):
    print("*** run end run(ERA)=", run.GetRunID())

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
  "My Event Action"

  #def BeginOfEventAction(self, event):
    #print("*** current event (BEA)=", event.GetEventID())
  #  pass

  #def EndOfEventAction(self, event):
    #print("*** current event (EEA)=", event.GetEventID())

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
  "My Stepping Action"

  def UserSteppingAction(self, step):
    #print("*** dE/dx in current step=", step.GetTotalEnergyDeposit())
    track= step.GetTrack()
    touchable= track.GetTouchable()
    pv= touchable.GetVolume()
    #print(pv.GetCopyNo())
    #print(touchable.GetReplicaNumber(0))

# ------------------------------------------------------------------
class MyField(G4MagneticField):
  "My Magnetic Field"

  def GetFieldValue(self, pos, time):
    bfield= G4ThreeVector()
    bfield.x= 0.
    bfield.y= 5.*tesla
    bfield.z= 0.
    return bfield
    
# ==================================================================
# main
# ==================================================================
qMaterials= gtest01.QMaterials()
qMaterials.Construct()

qDC= gtest01.QDetectorConstruction()
gRunManager.SetUserInitialization(qDC)

qPL= gtest01.QPhysicsList()
gRunManager.SetUserInitialization(qPL)

# set user actions...
#qPGA= gtest01.QPrimaryGeneratorAction()
myPGA= MyPrimaryGeneratorAction()
gRunManager.SetUserAction(myPGA)

#myRA= MyRunAction()
#gRunManager.SetUserAction(myRA)
  
myEA= MyEventAction()
gRunManager.SetUserAction(myEA)

mySA= MySteppingAction()
gRunManager.SetUserAction(mySA)

# set particle gun
#ApplyUICommand("/control/execute gun.mac")
#pg= qPGA.GetParticleGun()
pg= myPGA.particleGun
pg.SetParticleByName("e-")
pg.SetParticleEnergy(200.*MeV)
pg.SetParticlePosition(G4ThreeVector(0.,0.,-14.9)*cm)

# magnetic field
fieldMgr= gTransportationManager.GetFieldManager()
myField= G4UniformMagField(G4ThreeVector(0.,10.*tesla,0.))
#myField= MyField()
fieldMgr.SetDetectorField(myField)
fieldMgr.CreateChordFinder(myField)

gRunManager.Initialize()

# visualization
gControlExecute("vis.mac")
  
# beamOn
gRunManager.BeamOn(10)

