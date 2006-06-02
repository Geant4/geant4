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
#include "HistoManager.hh"
#include "EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
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
  (HistoManager::GetPointer())->BeginOfEvent();

  if(nSelected>0) {
    for(G4int i=0; i<nSelected; i++) {
      if(nEvt == selectedEvents[i]) {
        UI->ApplyCommand("/random/saveThisEvent");
      }
    }
  }

  // Switch on verbose mode
  //  if(hi->GetFirstEventToDebug() == nEvt) {
  if(9985 == nEvt) {
    verbose = 1;
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 2");

  }

  // Switch off verbose mode
  //  if(hi->GetLastEventToDebug() == nEvt-1) {
  if(nEvt == 9986) {
    verbose = 0;
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 0");
  }

  // Initialize user actions
  if(verbose > 0 || G4int(nEvt/printModulo)*printModulo == nEvt) {
    G4cout << "EventAction: Event # "
           << nEvt << " started" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt)
{
  (HistoManager::GetPointer())->EndOfEvent();
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    G4TrajectoryContainer* trjc = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trjc) n_trajectories = trjc->entries();

    for(G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* t = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") t->DrawTrajectory(1000);
      else if ((drawFlag == "charged")&&(t->GetCharge() != 0.))
                             t->DrawTrajectory(1000);
      else if ((drawFlag == "neutral")&&(t->GetCharge() == 0.))
                             t->DrawTrajectory(1000);
      else if ((drawFlag == "charged+n")&&((t->GetCharge() != 0.)||
                                           (t->GetCharge()==0.&&t->GetParticleName()=="neutron")))
                             t->DrawTrajectory(1000);
   }
  }

  if(verbose > 0) {
    G4cout << "EventAction: Event # "
           << nEvt << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
