//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALTBEventAction.cc,v 1.5 2003-02-14 15:55:20 pmendez Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE

#include "FCALTBEventAction.hh"

#include "FCALTBEventActionMessenger.hh"

#include "FCALRunAction.hh"

#include "FCALCalorHit.hh"

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

#include "FCALAnalysisManager.hh"

#include "G4ios.hh"
#include "iostream.h"
#include "fstream.h"

#include <AIDA/AIDA.h>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTBEventAction::FCALTBEventAction(FCALSteppingAction* SA)
  :calorimeterCollID(-1),drawFlag("all"),printModulo(10), StepAction(SA), eventMessenger(0)
{
  eventMessenger = new FCALTBEventActionMessenger(this);
  runManager = new FCALRunAction();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTBEventAction::~FCALTBEventAction()
{
  delete eventMessenger;
  eventMessenger = 0;
  delete  runManager;
  runManager = 0;
}

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
  G4cout << "N Tracks out of world " << NTracksOutOfWorld << G4endl;

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

  G4double OutOfWorld;

  OutTracks << NTracksOutOfWorld << G4endl;
  for(G4int i=1; i<= NTracksOutOfWorld ; i++){
    for(G4int j=1; j<11 ; j++) {
      //      G4double OutOfWorld = StepAction->GetOutOfWorldTracks(i,j);
      OutOfWorld = StepAction->GetOutOfWorldTracks(i,j);
      OutTracks << OutOfWorld << " " ; 
    }
    OutTracks << G4endl;

    G4double OutOfWorld2 = StepAction->GetOutOfWorldTracks(i,j);


#ifdef G4ANALYSIS_USE

    FCALAnalysisManager* analysis = FCALAnalysisManager::getInstance();
    analysis->analyseEnergyDep(OutOfWorld);

#endif

    //    RunAction->GetHisto(1)->fill(OutOfWorld);


//    RunAction->GetTuple(1)->fill(0,OutOfWorld);
//    RunAction->GetTuple(1)->fill(1,i);
//    RunAction->GetTuple(1)->fill(2,j);
//
//    
//    RunAction->GetTuple(1)->addRow();

  } 

  OutTracks.close();

  NSecondaries = StepAction->GetSecondaries(0,0);
  G4cout << "N Scondaries " << NSecondaries << G4endl;   
  
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
  
  G4double Secondary;
  
  SecndTracks << NSecondaries << G4endl;
  for(i=1; i<= NSecondaries ; i++){
    for(G4int j=1; j<11 ; j++) {
      Secondary = StepAction->GetSecondaries(i,j);
      SecndTracks << Secondary  << " " ; }
    SecndTracks << G4endl;
    G4double Secondary2 = StepAction->GetSecondaries(i,j);

#ifdef G4ANALYSIS_USE
    FCALAnalysisManager* analysis = FCALAnalysisManager::getInstance();
    analysis->analyseEnergyDep(Secondary);
#endif

//    RunAction->GetHisto(2)->fill(Secondary);
//
//    RunAction->GetTuple(2)->fill(0,Secondary);
//    RunAction->GetTuple(2)->fill(1,i);
//    RunAction->GetTuple(2)->fill(2,j);
//    
//    
//    RunAction->GetTuple(2)->addRow();
    

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
  EdepFCAL << G4endl;
  EdepFCAL.close();


#ifdef G4ANALYSIS_USE
  FCALAnalysisManager* analysis = FCALAnalysisManager::getInstance();
  analysis->analyseEnergyDep(EmEdep);
  analysis->analyseEnergyDep(HadEdep);

#endif

//  RunAction->GetHisto(3)->fill(EmEdep);
//  RunAction->GetHisto(4)->fill(HadEdep);
//
//
//  RunAction->GetTuple(3)->fill(0,EmEdep);
//  RunAction->GetTuple(3)->fill(1,HadEdep);
//  
//    
//    
//  RunAction->GetTuple(3)->addRow();
    


  G4cout << "Edep in FCAL1 FCAl2 : " << StepAction->GetEdepFCAL("FCALEm") << " ";
  G4cout << StepAction->GetEdepFCAL("FCALHad") << G4endl;
}

#endif
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
