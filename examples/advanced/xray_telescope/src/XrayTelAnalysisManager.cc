// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: XrayTelAnalysisManager.cc,v 1.1 2000-11-10 13:35:53 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Guy Barrand 14th Septembre 2000.

#ifdef G4ANALYSIS_USE

#include "g4std/fstream"

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4SteppingManager.hh"
/*
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

#include <assert.h>
#include "g4std/vector"
#include "g4std/fstream"
#include "g4std/iomanip"
#include "g4std/iostream"
#include "g4std/vector"
*/

#include <IHistogramFactory.h>

#ifdef G4ANALYSIS_CLOUD
#include <ICloudFactory.h>
#endif

#ifdef G4ANALYSIS_USE_LAB
#include "G4LabSystem.hh"
#endif
#ifdef G4ANALYSIS_USE_JAS
#include "G4JasSystem.hh"
#endif
#ifdef G4ANALYSIS_USE_LIZARD
#include "G4LizardSystem.hh"
#endif

#include "XrayTelAnalysisManager.hh"

XrayTelAnalysisManager::XrayTelAnalysisManager(
 const G4String& aSystem
)
:fDrawEvent(false)
,fEnteringEnergyHistogram(0)
,fYZ_Histogram(0)
#ifdef G4ANALYSIS_NON_AIDA
,fKineticCloud(0)
,fYZ_Cloud(0)
,fTuple(0)
#endif
{
#ifdef G4ANALYSIS_USE_LAB
  RegisterAnalysisSystem(new G4LabSystem);
#endif
#ifdef G4ANALYSIS_USE_JAS
  RegisterAnalysisSystem(new G4JasSystem);
#endif
#ifdef G4ANALYSIS_USE_LIZARD
  RegisterAnalysisSystem(new G4LizardSystem);
#endif
  //  The factory and histograms will be 
  // deleted by the analysis manager.
  IHistogramFactory* hfactory = GetHistogramFactory(aSystem);
  if(hfactory) {
    // Histogram creation :
    fEnteringEnergyHistogram = 
      hfactory->create1D("Entering energy",100,0,0.002);
    // Book the histogram for the 2D position. 
    // Instead of using a scatter plot, just book enough bins ...
    fYZ_Histogram = 
      hfactory->create2D("YZ",200, -50, 50, 200, -50, 50);
  }

  // Under study...
#ifdef G4ANALYSIS_NON_AIDA
  ICloudFactory* cloudFactory = GetCloudFactory(aSystem);
  if(cloudFactory) {
    fKineticCloud = cloudFactory->create1D("Kinetic Energy");
    fYZ_Cloud = cloudFactory->create2D("YZ");
  }
  fTuple = CreateTuple(aSystem,"XrayTel.root","XrayTel");
#endif
}
XrayTelAnalysisManager::~XrayTelAnalysisManager(){
#ifdef G4ANALYSIS_NON_AIDA
  delete fTuple;
#endif
}
void XrayTelAnalysisManager::BeginOfRun(){
  fEnteringEnergy.clear();
  fEnteringDirection.clear();

  if(fEnteringEnergyHistogram) 
    fEnteringEnergyHistogram->reset();
}
void XrayTelAnalysisManager::EndOfRun(){
  G4std::ofstream outscat("detector.hist", ios::app);
  G4cout << "End of Run summary" << G4endl << G4endl;
  G4double TotEnteringEnergy = 0.0;
  for (G4int i=0;i< fEnteringEnergy.size();i++)
    TotEnteringEnergy += fEnteringEnergy[i];
  G4cout << "Total Entering Detector : " << fEnteringEnergy.size()  << G4endl;
  G4cout << "Total Entering Detector Energy : " << TotEnteringEnergy  << G4endl;
  for (i=0;i<fEnteringEnergy.size();i++) {
    outscat << "  "
	    << fEnteringEnergy[i]
            << "  "
            << fEnteringDirection[i].x()
	    << "  "
            << fEnteringDirection[i].y()
	    << "  "
            << fEnteringDirection[i].z()
	    << G4endl;
  }
  outscat.close();

  Store();  
#ifdef G4ANALYSIS_NON_AIDA
  if(fTuple) fTuple->commit();
#endif
}
void XrayTelAnalysisManager::BeginOfEvent(){
  fDrawEvent = false;
}
void XrayTelAnalysisManager::EndOfEvent(const G4Event* aEvent, const G4String& aDrawFlag){
  if (fDrawEvent && aEvent){
    if ( G4VVisManager::GetConcreteInstance() ) {
      G4TrajectoryContainer* trajectoryContainer = aEvent->GetTrajectoryContainer();
      if ( trajectoryContainer ){ 
	G4int n_trajectories = trajectoryContainer->entries(); 
	for ( G4int i=0; i<n_trajectories; i++ ) {
	  G4Trajectory* trj = (G4Trajectory*)(*(aEvent->GetTrajectoryContainer()))[i];
	  if ( aDrawFlag == "all" ) trj->DrawTrajectory(50);
	  else if ( (aDrawFlag == "charged")&&(trj->GetCharge() > 0.) )
	    trj->DrawTrajectory(50); 
	  trj->ShowTrajectory(); 
	}
      }
    }
  }
}

void XrayTelAnalysisManager::Step(const G4SteppingManager* aSteppingManager) {
  if(!aSteppingManager) return;

  G4Track* track = aSteppingManager->GetTrack();
  G4Step* step = aSteppingManager->GetStep();
  G4int TrackID = track->GetTrackID();
  G4int StepNo = track->GetCurrentStepNumber();

  if(StepNo >= 10000) track->SetTrackStatus(fStopAndKill);

  G4String volName; 
  if ( track->GetVolume() ) 
    volName =  track->GetVolume()->GetName(); 
  G4String nextVolName;
  if ( track->GetNextVolume() ) 
    nextVolName =  track->GetNextVolume()->GetName();
 
  G4ThreeVector pos = track->GetPosition();

  //--- Entering Detector
  if(volName != "Detector_P" && nextVolName == "Detector_P") {

    fEnteringEnergy.push_back(track->GetKineticEnergy());
    fEnteringDirection.push_back(pos);

    if(fEnteringEnergyHistogram) {
      fEnteringEnergyHistogram->fill(track->GetKineticEnergy());
      Plot(fEnteringEnergyHistogram);
    }

    if(fYZ_Histogram) {
      fYZ_Histogram->fill(pos.y(), pos.z());
      Plot(fYZ_Histogram);
    }

#ifdef G4ANALYSIS_NON_AIDA
    if(fTuple) {
      fTuple->fill("EnteringEnergy",track->GetKineticEnergy());
      fTuple->fill("x",pos.x());
      fTuple->fill("y",pos.y());
      fTuple->fill("z",pos.z());
      fTuple->flush();
    }
    if(fKineticCloud) {
      fKineticCloud->fill(track->GetKineticEnergy());
      fAnalysisManager->Plot(fKineticCloud);
    }
    if(fYZ_Cloud) {
      fYZ_Cloud->fill(pos.y(),pos.z());
      fAnalysisManager->Plot(fYZ_Cloud);
    }
#endif

    fDrawEvent = true;
  }
} 

#endif
