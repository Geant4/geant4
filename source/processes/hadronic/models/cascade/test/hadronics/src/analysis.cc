#ifdef G4ANALYSIS_USE

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Step.hh"

#include "analysis.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"

analysis::analysis() {
  hfile  = new TFile("N02ROOT.root","RECREATE","Demo");
  histoTrajectories = new TH1F("histoTrajectories" , "Number of trajectories / event", 30 , 0, 30);
  histoE = new TH1F("histoE" , "Energy deposited / step", 100, 0 , 1);
  ntuple = new TNtuple("ntuple" , "Demo ntuple","runId:eventId:nTrajectories:random:i");
}

analysis::~analysis() {
  hfile->Close();

  delete histoTrajectories;
  delete histoE;
  delete ntuple;
  delete hfile;
}

void analysis::BeginOfRun(const G4Run* aRun){
  runId = aRun->GetRunID();
  G4cout << "### Run " << runId << " start." << G4endl;
}

void analysis::EndOfRun(const G4Run*){
  hfile->Write(); 
}

void analysis::BeginOfEvent(const G4Event*){
}

void analysis::EndOfEvent(const G4Event* anEvent){

  G4int eventId = anEvent->GetEventID();
 
  G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
  G4int nTrajectories = 0;
  if (trajectoryContainer) nTrajectories = trajectoryContainer->entries();

  G4PrimaryParticle* primary = anEvent->GetPrimaryVertex(0)->GetPrimary(0);;
  histoTrajectories->Fill(nTrajectories);
  ntuple->Fill(runId, eventId, nTrajectories, 0, 0);
}

void analysis::Step(const G4Step* aStep){
  G4double edep = aStep->GetTotalEnergyDeposit()/ MeV;
  histoE->Fill(edep);
}

#endif
