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
// $Id: Tst33VisEventAction.cc,v 1.3 2002-11-04 10:57:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst33VisEventAction.hh"

#include "Tst33VisEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst33VisEventAction::Tst33VisEventAction()
:drawFlag("all"),printModulo(1),
 eventMessenger(0)
{
  eventMessenger = new Tst33VisEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst33VisEventAction::~Tst33VisEventAction()
{
  delete eventMessenger;
}

void Tst33VisEventAction::Clear() {
  
}


void Tst33VisEventAction::SetCell_19_Scorer(const G4CellScorer *scorer){
  G4cout << "Tst33VisEventAction::SetCell_19_Scorer: no action taken!" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst33VisEventAction::BeginOfEventAction(const G4Event* evt)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst33VisEventAction::EndOfEventAction(const G4Event* evt)
{
  
  if (G4VVisManager::GetConcreteInstance()) {
    G4TrajectoryContainer* trajectoryContainer(0);
    trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) {
      n_trajectories = trajectoryContainer->entries();
    }

    for (G4int i=0; i<n_trajectories; i++)  { 
      G4Trajectory* trj = 0;
      trj = dynamic_cast<G4Trajectory*>((*trajectoryContainer)[i]);
      if (trj) {
	G4bool charged(false);
	charged = G4std::fabs(trj->GetCharge()) > 0;
	if (drawFlag == "all") {
	  trj->DrawTrajectory(50);
	}
	else if ((drawFlag == "charged")&&(charged)) {
	  trj->DrawTrajectory(50);
	}
	else if ((drawFlag == "neutral")&&(!charged)) {
	  trj->DrawTrajectory(50);
	}
      }
      else {
	G4cerr << " Tst33VisEventAction::EndOfEventAction: failed to dynamic cast to G4Trajectory!" << G4endl;
      }
    }
  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
