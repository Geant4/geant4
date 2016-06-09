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
// $Id: EventAction.cc,v 1.6 2010-06-07 05:40:46 perl Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "EventAction.hh"
#include "HistoManager.hh"
#include "EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4Gamma.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction():
  nEvt(0),
  printModulo(100),
  verbose(0),
  drawFlag("all")
{
  eventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*)
{
  // New event
  nEvt++;
  HistoManager* hi = HistoManager::GetPointer();
  hi->BeginOfEvent();

  // Switch on verbose mode

  if(hi->FirstEventToDebug() == nEvt) {
    verbose = 2;
    hi->SetVerbose(2);
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 2");
  }

  // Switch off verbose mode
  if(hi->LastEventToDebug() == nEvt-1) {
    verbose = 0;
    hi->SetVerbose(0);
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 0");
  }


  // Initialize user actions
  if(verbose > 0) {
    G4cout << "EventAction: Event # "
           << nEvt << " started" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event*)
{
  (HistoManager::GetPointer())->EndOfEvent();

  if(verbose > 0) {
    G4cout << "EventAction: Event # "
           << nEvt << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
