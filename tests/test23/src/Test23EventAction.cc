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
//
// $Id: Test23EventAction.cc,v 1.1 2004-03-18 11:02:26 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- Test23EventAction class ----------------
//                 by Mikhail Kossov, December 2003.
//  Test23EventAction class of the CHIPS Simulation Branch in GEANT4

#include "Test23EventAction.hh"
#include "Test23EventActionMessenger.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4Gamma.hh"
#include "G4ios.hh"


Test23EventAction::Test23EventAction():
  nEvt(0),
  printModulo(100),
  verbose(0),
  drawFlag("all")
{
  eventMessenger = new Test23EventActionMessenger(this);
}

Test23EventAction::~Test23EventAction()
{}

void Test23EventAction::BeginOfEventAction(const G4Event*)
{
  // New event
  nEvt++;
  if(100*(nEvt/100)==nEvt) G4cout<<"Test23 event # "<<nEvt<<G4endl;

  // Switch on verbose mode
  // Initialize user actions
  if(verbose > 0) {
    G4cout << "EventAction: Event # "
           << nEvt << " started" << G4endl;
  }

}

//void Test23EventAction::EndOfEventAction(const G4Event* evt)
void Test23EventAction::EndOfEventAction(const G4Event* )
{
  //(HistoManager::GetPointer())->EndOfEvent();
  //G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  //if(pVVisManager) {
  //  G4TrajectoryContainer* trjc = evt->GetTrajectoryContainer();
  //  G4int n_trajectories = 0;
  //  if (trjc) n_trajectories = trjc->entries();

  //  for(G4int i=0; i<n_trajectories; i++) {
  //    G4Trajectory* t = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
  //    if (drawFlag == "all") t->DrawTrajectory(1000);
  //    else if ((drawFlag == "charged")&&(t->GetCharge() != 0.))
  //                           t->DrawTrajectory(1000);
  //  }
  //}

  if(verbose > 0) G4cout << "EventAction: Event # " << nEvt << " ended" << G4endl;
}
