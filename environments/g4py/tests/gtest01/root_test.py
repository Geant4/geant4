#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest01
#   - check basic control flow
# ==================================================================
from Geant4 import *
import gtest01
import ROOT

# ==================================================================
#                          ROOT PART                               #
# ==================================================================

# ------------------------------------------------------------------
def init_root():
# ------------------------------------------------------------------  
  ROOT.gROOT.Reset()

  # plot style
  ROOT.gStyle.SetTextFont(82)
  ROOT.gStyle.SetTitleFont(82, "X")
  ROOT.gStyle.SetLabelFont(82, "X")
  ROOT.gStyle.SetTitleFont(82, "Y")
  ROOT.gStyle.SetLabelFont(82, "Y")

  #ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetErrorX(0)

  canvas= ROOT.TCanvas("g4py_plots",
                       "Geant4Py Sample Plots",
                       620, 30, 600, 400)

  #canvas.Divide(2,2);
  #canvas.SetFillColor(29)

  canvas.SetGrid()

  return canvas

# ------------------------------------------------------------------
def hini():
# ------------------------------------------------------------------  
  global hist1
  hist1= ROOT.TH1D("dE/dx/step", "dE/dx", 100, 0., 2000.)
  hist1.SetXTitle("(keV)")


# ------------------------------------------------------------------
def hshow():
# ------------------------------------------------------------------  
  hist1.Draw()
  
# ==================================================================
#                         Geant4 PART                              #
# ==================================================================

# ==================================================================
#   user actions in python
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

  def BeginOfRunAction(self, run):
    print "*** #event to be processed (BRA)=",
    run.numberOfEventToBeProcessed

  def EndOfRunAction(self, run):
    print "*** run end run(ERA)=", run.runID

# ------------------------------------------------------------------
class MyEventAction(G4UserEventAction):
  "My Event Action"

  def BeginOfEventAction(self, event):
    print "*** current event (BEA)=", event.eventID

  def EndOfEventAction(self, event):
    print "*** current event (EEA)=", event.eventID

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
  "My Stepping Action"

  def UserSteppingAction(self, step):
    #print "*** dE/dx in current step=", step.GetTotalEnergyDeposit()
    dedx= step.GetTotalEnergyDeposit()
    if(dedx>0):
      hist1.Fill(dedx/HEPUnit.keV)

# ==================================================================
#   main
# ==================================================================
g4pyCanvas= init_root()
hini()

app= gtest01.MyApplication()
app.Configure()

# set user actions...
myPGA= MyPrimaryGeneratorAction()
gRunManager.SetUserAction(myPGA)

myRA= MyRunAction()
gRunManager.SetUserAction(myRA)
  
#myEA= MyEventAction()
#gRunManager.SetUserAction(myEA)

mySA= MySteppingAction()
gRunManager.SetUserAction(mySA)

  
# set particle gun
#ApplyUICommand("/control/execute gun.mac")
pg= myPGA.particleGun
pg.SetParticleByName("e-")
pg.SetParticleEnergy(200.*HEPUnit.MeV)
pg.SetParticleMomentumDirection(G4ThreeVector(0.2, 0., 1.))
pg.SetParticlePosition(G4ThreeVector(0.,0.,-14.9)*HEPUnit.cm)
  
# visualization
ApplyUICommand("/control/execute vis.mac")
  
# beamOn
gRunManager.BeamOn(1000)

#
hshow()


