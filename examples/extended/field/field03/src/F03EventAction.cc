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
/// \file field/field03/src/F03EventAction.cc
/// \brief Implementation of the F03EventAction class
//
//
// $Id: F03EventAction.cc 76602 2013-11-13 08:33:35Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F03EventAction.hh"

#include "F03RunAction.hh"

#include "F03CalorHit.hh"
#include "F03EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03EventAction::F03EventAction(F03RunAction* action)
 : G4UserEventAction(),
   fCalorimeterCollID(-1),
   fEventMessenger(0),
   fRunAction(action),
   fVerboseLevel(0),
   fPrintModulo(10000)
{
  fEventMessenger = new F03EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F03EventAction::~F03EventAction()
{
  delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();
  if (evtNb%fPrintModulo == 0)
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;

  if (fVerboseLevel>1)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;

  if (fCalorimeterCollID==-1)
    {
     G4SDManager * sdManager = G4SDManager::GetSDMpointer();
     fCalorimeterCollID = sdManager->GetCollectionID("CalCollection");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F03EventAction::EndOfEventAction(const G4Event* evt)
{
  if (fVerboseLevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;

  // save rndm status
  if (fRunAction->GetRndmFreq() == 2)
    {
     G4Random::saveEngineStatus("endOfEvent.rndm");
     G4int evtNb = evt->GetEventID();
     if (evtNb%fPrintModulo == 0)
       {
         G4cout << "\n---> End of Event: " << evtNb << G4endl;
         G4Random::showEngineStatus();
       }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
