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
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17EventAction.hh"

#include "Test17RunAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17EventAction::Test17EventAction(Test17RunAction* Test17RA)
 :runaction(Test17RA),
  evtNo(-1),
  good(false)
{
  verbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17EventAction::~Test17EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::BeginOfEventAction(const G4Event* evt)
{
  evtNo = evt->GetEventID();
  edep = 0.0;
  nCharged = 0;
  nNeutral = 0;
  totLength = 0.0;

  if(verbose>0)
    G4cout << "<<< Event  " << evtNo << " started." << G4endl;

  /*
  if(14 == evtNo) {
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 2");
    (G4UImanager::GetUIpointer())->ApplyCommand("/stepping/verbose 2");
  }
  if(16 == evtNo) {
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 0");
    (G4UImanager::GetUIpointer())->ApplyCommand("/stepping/verbose 0");
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventAction::EndOfEventAction(const G4Event*)
{
   // count event, add deposits to the sum ...
  if(good) {
    runaction->AddTrackLength(totLength) ;
    runaction->AddEdep(edep);
    runaction->CountParticles(nCharged,nNeutral);
    runaction->CountEvent();
  }
  if(verbose>0) {
      G4cout << " Ncharged= " << nCharged 
             << "; Nneutral= " << nNeutral
             << "; Length(mm)= " << totLength/mm 
             << "; Edep(MeV)= " << edep/MeV 
             << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  















