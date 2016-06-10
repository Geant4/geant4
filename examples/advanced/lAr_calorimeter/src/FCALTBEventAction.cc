//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: FCALTBEventAction.cc 73002 2013-08-15 08:09:37Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FCALTBEventAction.hh"
#include "FCALTBEventActionMessenger.hh"
#include "FCALRunAction.hh"
#include "FCALCalorHit.hh"
#include <vector>

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "FCALSteppingAction.hh"

#include "FCALAnalysisManager.hh"

#include "G4ios.hh"
#include <fstream>
#include <iostream>

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTBEventAction::FCALTBEventAction(FCALSteppingAction* SA)
  :drawFlag("all"),printModulo(10), StepAction(SA), eventMessenger(0)
{
  eventMessenger = new FCALTBEventActionMessenger(this);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTBEventAction::~FCALTBEventAction()
{
  delete eventMessenger;
  eventMessenger = 0;
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

void FCALTBEventAction::EndOfEventAction(const G4Event*)
{
  // Fill histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  G4int i,j;
  NTracksOutOfWorld = StepAction->GetOutOfWorldTracks(0, 0); 
  G4cout << "N Tracks out of world " << NTracksOutOfWorld << G4endl;

  // Write Leaving Particles in File
  //--------------------------------
  G4String FileName1 = "OutTracks_802_1mm.dat";
  std::ios::openmode iostemp1;
  if(Init1 == 0) {
    iostemp1 = std::ios::out;
    Init1++;
  } else {
    iostemp1 = std::ios::out|std::ios::app; // std::ios::app;  
  };
  std::ofstream OutTracks(FileName1, iostemp1);

  OutTracks << NTracksOutOfWorld << G4endl;

  G4double OutOfWorld;
  for(i=1; i<= NTracksOutOfWorld ; i++){
    for(j=1; j<11 ; j++) {
      //      G4double OutOfWorld = StepAction->GetOutOfWorldTracks(i,j);
      OutOfWorld = StepAction->GetOutOfWorldTracks(i,j);
      OutTracks << OutOfWorld << " " ; 
    }
    OutTracks << std::endl;
    // G4double OutOfWorld2 = StepAction->GetOutOfWorldTracks(i,j);
  } 
  OutTracks.close();

  man->FillH1(1,NTracksOutOfWorld);

  NSecondaries = StepAction->GetSecondaries(0,0);
  G4cout << "N Scondaries " << NSecondaries << G4endl;   
  
  // Write Secondary Particles in File
  //----------------------------------
  G4String FileName2 = "SecndTracks_802_1mm.dat";
  std::ios::openmode iostemp2;
  if(Init2 == 0) {
    iostemp2 = std::ios::out;
    Init2++;
  } else {
    iostemp2 = std::ios::out|std::ios::app; // std::ios::app;  
  };
  
  std::ofstream SecndTracks(FileName2, iostemp2);
  
  SecndTracks << NSecondaries << std::endl;

  G4double Secondary;  
  for(i=1; i<= NSecondaries ; i++){
    for(j=1; j<11 ; j++) {
      Secondary = StepAction->GetSecondaries(i,j);
      SecndTracks << Secondary  << " " ; 
    }
    SecndTracks << std::endl;
    // G4double Secondary2 = StepAction->GetSecondaries(i,j);
  }
  SecndTracks.close();
  
  man->FillH1(2,NSecondaries);

  // Write Edep in FCAL1 and FCAL2 
  G4String FileName3 = "EdepFCAL_802_1mm.dat";
  std::ios::openmode iostemp3;
  if(Init3 == 0) {
    iostemp3 = std::ios::out;
    Init3++;
  } else {
    iostemp3 = std::ios::out|std::ios::app; // std::ios::app;  
  };
  
  std::ofstream EdepFCAL(FileName3, iostemp3);
  
  G4double EmEdep  = StepAction->GetEdepFCAL("FCALEm");
  G4double HadEdep = StepAction->GetEdepFCAL("FCALHad");

  EdepFCAL << EmEdep << " ";
  EdepFCAL << HadEdep; 
  EdepFCAL << std::endl;
  EdepFCAL.close();

  man->FillH1(3,EmEdep/MeV);
  man->FillH1(4,HadEdep/MeV);


  G4cout<<"EmEdep is="<<EmEdep/MeV << " MeV " << G4endl;
  G4cout<<"HadEdep is="<<HadEdep << " MeV" << G4endl;

  G4cout << "Edep in FCAL1 FCAl2 : " << StepAction->GetEdepFCAL("FCALEm") << " ";
  G4cout << StepAction->GetEdepFCAL("FCALHad") << G4endl;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
