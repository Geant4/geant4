// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALTBEventAction.cc,v 1.2 2002-10-03 11:34:42 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALTBEventAction.hh"

#include "FCALCalorHit.hh"

#ifdef G4ANALYSIS_USE
#include "FCALAnalysisManager.hh"
#endif

//#include "g4rw/tvordvec.h"
#include "g4std/vector"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "FCALSteppingAction.hh"

#include "G4ios.hh"
#include "iostream.h"
#include "fstream.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTBEventAction::FCALTBEventAction(FCALSteppingAction* SA)
:calorimeterCollID(-1),drawFlag("all"),printModulo(10), StepAction(SA)
{;}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTBEventAction::~FCALTBEventAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTBEventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
    if (evtNb%printModulo == 0)
    { 
      G4cout << "\n---> Begin of event: " << evtNb+1 << G4endl;
//      HepRandom::showEngineStatus();
    }
    
  NTracksOutOfWorld = 0;
  NSecondaries = 0;
 
  StepAction->initialize(evtNb+1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTBEventAction::EndOfEventAction(const G4Event* evt)
{
  NTracksOutOfWorld = StepAction->GetOutOfWorldTracks(0, 0); 
  G4cout << "N Tracks out of world " << NTracksOutOfWorld << endl;

  // Write Leaving Particles in File
  //--------------------------------
  G4String FileName1 = "OutTracks_802_1mm.dat";
  G4int iostemp1;
  if(Init1 == 0) {
    iostemp1 = ios::out;
    Init1++;
  } else {
    iostemp1 = ios::out|ios::app; // ios::app;  
  };
  ofstream OutTracks(FileName1, iostemp1);

  OutTracks << NTracksOutOfWorld << endl;
  for(G4int i=1; i<= NTracksOutOfWorld ; i++){
    for(G4int j=1; j<11 ; j++) {
      G4double OutOfWorld = StepAction->GetOutOfWorldTracks(i,j);
      OutTracks << OutOfWorld << " " ; }
    OutTracks << endl;
#ifdef G4ANALYSIS_USE
	// pass first event for histogram record:
      G4double OutOfWorld2 = StepAction->GetOutOfWorldTracks(i,j);
	FCALAnalysisManager* analysis = FCALAnalysisManager::getInstance();
	analysis->NumOutOfWorld(OutOfWorld2,i,j);
#endif
  }
  OutTracks.close();

  NSecondaries = StepAction->GetSecondaries(0,0);
  G4cout << "N Scondaries " << NSecondaries << endl;   



  
  // Write Secondary Particles in File
  //--------------------------------
  G4String FileName2 = "SecndTracks_802_1mm.dat";
  G4int iostemp2;
  if(Init2 == 0) {
    iostemp2 = ios::out;
    Init2++;
  } else {
    iostemp2 = ios::out|ios::app; // ios::app;  
  };

  ofstream SecndTracks(FileName2, iostemp2);

  SecndTracks << NSecondaries << endl;
  for(i=1; i<= NSecondaries ; i++){
    for(G4int j=1; j<11 ; j++) {
      G4double Secondary = StepAction->GetSecondaries(i,j);
      SecndTracks << Secondary  << " " ; }
    SecndTracks << endl;
#ifdef G4ANALYSIS_USE
	// pass first event for histogram record:
      G4double Secondary2 = StepAction->GetSecondaries(i,j);
	FCALAnalysisManager* analysis = FCALAnalysisManager::getInstance();
	analysis->Secondaries(Secondary2,i,j);
#endif
  }
  SecndTracks.close();


  // Write Edep in FCAL1 and FCAL2 
  G4String FileName3 = "EdepFCAL_802_1mm.dat";
  G4int iostemp3;
  if(Init3 == 0) {
    iostemp3 = ios::out;
    Init3++;
  } else {
    iostemp3 = ios::out|ios::app; // ios::app;  
  };

  ofstream EdepFCAL(FileName3, iostemp3);
  
  G4double EmEdep  = StepAction->GetEdepFCAL("FCALEm");
  G4double HadEdep = StepAction->GetEdepFCAL("FCALHad");

  EdepFCAL << EmEdep << " ";
  EdepFCAL << HadEdep; 
  EdepFCAL << endl;
  EdepFCAL.close();

#ifdef G4ANALYSIS_USE
	// pass first event for histogram record:
	FCALAnalysisManager* analysis = FCALAnalysisManager::getInstance();
	analysis->Edep(EmEdep,HadEdep);
#endif



  G4cout << "Edep in FCAL1 FCAl2 : " << StepAction->GetEdepFCAL("FCALEm") << " ";
  G4cout << StepAction->GetEdepFCAL("FCALHad") << endl;
}


  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
