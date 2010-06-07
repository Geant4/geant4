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
// -------------------------------------------------------------
//
//
//      ---------- EventAction -------------
//
//  Modified: 05.04.01 Vladimir Ivanchenko new design of IBREM
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "EventAction.hh"
#include "Histo.hh"
#include "EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction():
  nEvt(0),
  verbose(0),
  printModulo(100),
  nSelected(0),
  //  drawFlag("all")
  drawFlag("charged")
  // drawFlag("neutral")
{
  // drawFlags = all, charged, neutral, charged+n
  eventMessenger = new EventActionMessenger(this);
  UI = G4UImanager::GetUIpointer();
  selectedEvents.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*)
{
  // New event
  nEvt++;
  (Histo::GetPointer())->AddEvent();

  if(nSelected>0) {
    for(G4int i=0; i<nSelected; i++) {
      if(nEvt == selectedEvents[i]) {
        UI->ApplyCommand("/random/saveThisEvent");
      }
    }
  }

  // Switch on verbose mode
  /*
  if(hi->GetFirstEventToDebug() == nEvt) {
    verbose = 2;
    hi->SetVerbose(2);
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 2");

    const G4ProcessManager* pm = G4Gamma::Gamma()->GetProcessManager();
    const G4ProcessVector* pv = pm->GetProcessList();
    G4int np = pm->GetProcessListLength();
    for(G4int i=0; i<np; i++) {
     ((*pv)[i])->SetVerboseLevel(2);
    }

  }

  // Switch off verbose mode
  if(hi->GetLastEventToDebug() == nEvt-1) {
    verbose = 0;
    hi->SetVerbose(0);
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 0");
  }
  */

  // Initialize user actions
  if(verbose > 0) {
    G4cout << "EventAction: Event # "
           << nEvt << " started" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event*)
{
  (Histo::GetPointer())->SaveEvent();

  if(verbose > 0) {
    G4cout << "EventAction: Event # "
           << nEvt << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
