// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01AnalysisManager.cc,v 1.7 2000-11-10 14:45:33 gbarrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Guy Barrand 14th Septembre 2000.

#ifdef G4ANALYSIS_USE

#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"

#include <IHistogramFactory.h>

#ifdef G4ANALYSIS_USE_LAB
#include "G4LabSystem.hh"
#endif
#ifdef G4ANALYSIS_USE_JAS
#include "G4JasSystem.hh"
#endif

#include "AnaEx01CalorHit.hh"
#include "AnaEx01AnalysisManager.hh"

AnaEx01AnalysisManager::AnaEx01AnalysisManager(
  const G4String& aSystem
)
:fCalorimeterCollID(-1)
,fEAbs(0)
,fLAbs(0)
,fEGap(0)
,fLGap(0)
{
#ifdef G4ANALYSIS_USE_LAB
  RegisterAnalysisSystem(new G4LabSystem);
#endif
#ifdef G4ANALYSIS_USE_JAS
  RegisterAnalysisSystem(new G4JasSystem);
#endif
  // The factory and histograms will be deleted by the analysis manager.
  IHistogramFactory* hfactory = GetHistogramFactory(aSystem);
  if(!hfactory) return;
  fEAbs = hfactory->create1D("EAbs",100,0,100);
  fLAbs = hfactory->create1D("LAbs",100,0,100);
  fEGap = hfactory->create1D("EGap",100,0,10);
  fLGap = hfactory->create1D("LGap",100,0,100);
}
void AnaEx01AnalysisManager::BeginOfRun(const G4Run* aRun){
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void AnaEx01AnalysisManager::EndOfRun(const G4Run*){
  Store();
  //if(fEAbs) Plot(fEAbs);
}

void AnaEx01AnalysisManager::BeginOfEvent(const G4Event*){
  if(fCalorimeterCollID==-1) {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    fCalorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 
}

void AnaEx01AnalysisManager::EndOfEvent(const G4Event* aEvent){
  if(!fEAbs) return; // No histo booked !
  //if(!fTuple) return; // No tuple booked !

  G4int evtNb = aEvent->GetEventID();

  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  AnaEx01CalorHitsCollection* CHC = 
    HCE ? (AnaEx01CalorHitsCollection*)(HCE->GetHC(fCalorimeterCollID)) : 0;

  if(CHC) {
    G4int n_hit = CHC->entries();
    //G4double totEAbs=0, totLAbs=0, totEGap=0, totLGap=0;
    for (G4int i=0;i<n_hit;i++) {
      //totEAbs += (*CHC)[i]->GetEdepAbs(); 
      //totLAbs += (*CHC)[i]->GetTrakAbs();
      //totEGap += (*CHC)[i]->GetEdepGap(); 
      //totLGap += (*CHC)[i]->GetTrakGap();
      fEAbs->fill((*CHC)[i]->GetEdepAbs());
      fLAbs->fill((*CHC)[i]->GetTrakAbs());
      fEGap->fill((*CHC)[i]->GetEdepGap());
      fLGap->fill((*CHC)[i]->GetTrakGap());
    }
  }	
  
}
void AnaEx01AnalysisManager::Step(const G4Step*){}

#endif
