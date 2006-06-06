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
// $Id: EventAction.cc,v 1.2 2006-06-06 19:48:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// EventAction
//
// Created: 31.04.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
// 

#include "EventAction.hh"
#include "G4Event.hh"
#include "HistoManager.hh"
#include "EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction():
  printModulo(100),
  nSelected(0),
  drawFlag("all"),
  debugStarted(false)
{
  eventMessenger = new EventActionMessenger(this);
  UI = G4UImanager::GetUIpointer();
  selectedEvents.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  // New event
  G4int nEvt = evt->GetEventID();

  if(nSelected>0) {
    for(G4int i=0; i<nSelected; i++) {
      if(nEvt == selectedEvents[i]) {
        UI->ApplyCommand("/random/saveThisEvent");
        UI->ApplyCommand("/tracking/verbose  2");
        debugStarted = true;
        break;
      }
    }
  }

  // Initialize user actions
  if(HistoManager::GetPointer()->GetVerbose() > 0 || 
     G4int(nEvt/printModulo)*printModulo == nEvt) 
    G4cout << "EventAction: Event # "
           << nEvt << " started" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt)
{
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

  if(debugStarted) {
    UI->ApplyCommand("/tracking/verbose  0");
    debugStarted = false;
  }

  if(HistoManager::GetPointer()->GetVerbose() > 1) 
    G4cout << "EventAction: Event ended" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
