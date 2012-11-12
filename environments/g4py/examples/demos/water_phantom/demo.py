#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#   gtest01
#   - check basic control flow
# ==================================================================
from Geant4 import *
import demo_wp
import g4py.MedicalBeam
import ROOT

# ==================================================================
#                          ROOT PART                               #
# ==================================================================

# ------------------------------------------------------------------
def init_root():
# ------------------------------------------------------------------  
  ROOT.gROOT.Reset()

  # plot style
  ROOT.gStyle.SetTextFont(42)
  ROOT.gStyle.SetTitleFont(42, "X")
  ROOT.gStyle.SetLabelFont(42, "X")
  ROOT.gStyle.SetTitleFont(42, "Y")
  ROOT.gStyle.SetLabelFont(42, "Y")
  
  global gCanvas
  gCanvas= ROOT.TCanvas("water_phantom_plots",
                        "Water Phantom Demo Plots",
                        620, 30, 800, 800)

# ------------------------------------------------------------------
def hini():
# ------------------------------------------------------------------  
  global gPad1
  gPad1= ROOT.TPad("2D", "2D", 0.02, 0.5, 0.98, 1.)
  gPad1.Draw()
  gPad1.cd()

  ROOT.gStyle.SetPalette(1);
   
  global hist_dose2d
  hist_dose2d= ROOT.TH2D("2D Dose", "Dose Distribution",
                         200, 0., 400.,
                         81, -81., 81.)
  hist_dose2d.SetXTitle("Z (mm)")
  hist_dose2d.SetYTitle("X (mm)")
  hist_dose2d.SetStats(0)
  hist_dose2d.Draw("colz")
  
  gCanvas.cd()
  global gPad2
  gPad2= ROOT.TPad("Z", "Z", 0.02, 0., 0.98, 0.5)
  gPad2.Draw()
  gPad2.cd()
  
  global hist_dosez
  hist_dosez= ROOT.TH1D("Z Dose", "Depth Dose", 200, 0., 400.)
  hist_dosez.SetXTitle("(mm)")
  hist_dosez.SetYTitle("Accumulated Dose (MeV)")
  hist_dosez.Draw()

# ------------------------------------------------------------------
def hshow():
# ------------------------------------------------------------------
  gPad1.cd()
  hist_dose2d.Draw("colz")
  gPad2.cd()
  hist_dosez.Draw()
  
  
# ==================================================================
#                         Geant4 PART                              #
# ==================================================================

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
    gPad1.Modified()
    gPad1.Update()
    gPad2.Modified()
    gPad2.Update()
    ROOT.gSystem.ProcessEvents()

# ------------------------------------------------------------------
class MySteppingAction(G4UserSteppingAction):
  "My Stepping Action"

  def UserSteppingAction(self, step):
    pass
  
# ------------------------------------------------------------------
class ScoreSD(G4VSensitiveDetector):
  "SD for score voxels"

  def __init__(self):
    G4VSensitiveDetector.__init__(self, "ScoreVoxel")

  def ProcessHits(self, step, rohist):
    preStepPoint= step.GetPreStepPoint()
    if(preStepPoint.GetCharge() == 0):
      return

    track= step.GetTrack()
    touchable= track.GetTouchable()
    voxel_id= touchable.GetReplicaNumber()
    dedx= step.GetTotalEnergyDeposit()
    xz= posXZ(voxel_id)
    hist_dose2d.Fill(xz[1], xz[0], dedx/MeV)
    if( abs(xz[0]) <= 100 ):
      hist_dosez.Fill(xz[1],  dedx/MeV)

# ------------------------------------------------------------------
def posXZ(copyN):
  dd= 2.*mm
  nx= 81

  iz= copyN/nx
  ix= copyN-iz*nx-nx/2
  
  x0= ix*dd
  z0= (iz+0.5)*dd
  return (x0,z0)


# ==================================================================
#   main
# ==================================================================
# init ROOT...
init_root()
hini()

# configure application
#app= demo_wp.MyApplication()
#app.Configure()

myMaterials= demo_wp.MyMaterials()
myMaterials.Construct()

myDC= demo_wp.MyDetectorConstruction()
gRunManager.SetUserInitialization(myDC)

myPL= FTFP_BERT()
gRunManager.SetUserInitialization(myPL)

# set user actions...
myPGA= MyPrimaryGeneratorAction()
gRunManager.SetUserAction(myPGA)

myRA= MyRunAction()
gRunManager.SetUserAction(myRA)
  
myEA= MyEventAction()
gRunManager.SetUserAction(myEA)

#mySA= MySteppingAction()
#gRunManager.SetUserAction(mySA)

# set particle gun
#pg= myPGA.particleGun
#pg.SetParticleByName("proton")
#pg.SetParticleEnergy(230.*MeV)
#pg.SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.))
#pg.SetParticlePosition(G4ThreeVector(0.,0.,-50.)*cm)

# medical beam
beam= g4py.MedicalBeam.Construct()
beam.particle= "proton"
beam.kineticE= 230.*MeV
#beam.particle= "gamma"
#beam.kineticE= 1.77*MeV
beam.sourcePosition= G4ThreeVector(0.,0.,-100.*cm)
beam.SSD= 100.*cm
beam.fieldXY= [5.*cm, 5.*cm]
  
# initialize
gRunManager.Initialize()

# set SD (A SD should be set after geometry construction)
scoreSD= ScoreSD()
myDC.SetSDtoScoreVoxel(scoreSD)

# visualization
gApplyUICommand("/control/execute vis.mac")
  
# beamOn
gRunManager.BeamOn(100)

#ROOT.gSystem.Run()

