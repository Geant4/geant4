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
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestEventAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestEventAction.hh"
#include "hTestHisto.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::hTestEventAction(const hTestDetectorConstruction* det):
  theDet(det),
  nEvt(0),
  verbose(0),
  drawFlag("all")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::~hTestEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::BeginOfEventAction(const G4Event* evt)
{  
  // New event
  if(nEvt == 0) (hTestHisto::GetPointer())->BeginOfHisto(0);
  nEvt++;

  // Switch on verbose mode
  if(theDet->GetFirstEventToDebug() == nEvt) {
    verbose = 2;
    (hTestHisto::GetPointer())->SetVerbose(2);    
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 2");
    (G4UImanager::GetUIpointer())->ApplyCommand("/stepping/verbose 2");
  }
    
  // Switch off verbose mode
  if(theDet->GetLastEventToDebug() == nEvt-1) {
    verbose = 0;
    (hTestHisto::GetPointer())->SetVerbose(0);    
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 0");
    (G4UImanager::GetUIpointer())->ApplyCommand("/stepping/verbose 0");
  }

  // Initialize user actions
  if(verbose > 0) {
    G4cout << "hTestEventAction: Event # " 
           << nEvt << " started" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::EndOfEventAction(const G4Event* evt)
{

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    G4TrajectoryContainer* trjc = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trjc) n_trajectories = trjc->entries();  

    for(G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* t = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") t->DrawTrajectory(50);
      else if ((drawFlag == "charged")&&(t->GetCharge() != 0.))
                             t->DrawTrajectory(50); 
    }
  }  

  if(verbose > 0) {
    G4cout << "hTestEventAction: Event # " 
           << nEvt << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  













