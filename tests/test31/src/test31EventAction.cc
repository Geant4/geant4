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
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31EventAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31EventAction.hh"
#include "test31Histo.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31EventAction::test31EventAction(const test31DetectorConstruction* det):
  theDet(det),
  nEvt(0),
  verbose(0),
  drawFlag("all")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31EventAction::~test31EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31EventAction::BeginOfEventAction(const G4Event* evt)
{  
  // New event
  nEvt++;
  (test31Histo::GetPointer())->AddEvent();

  // Switch on verbose mode
  if(theDet->GetFirstEventToDebug() == nEvt) {
    verbose = 2;
    (test31Histo::GetPointer())->SetVerbose(2);    
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 2");
    /*
    const G4ProcessManager* pm = G4Gamma::Gamma()->GetProcessManager();
    const G4ProcessVector* pv = pm->GetProcessList();
    G4int np = pm->GetProcessListLength();
    for(G4int i=0; i<np; i++) {
     ((*pv)[i])->SetVerboseLevel(2);
    }
    */
  }
    
  // Switch off verbose mode
  if(theDet->GetLastEventToDebug() == nEvt-1) {
    verbose = 0;
    (test31Histo::GetPointer())->SetVerbose(0);    
    (G4UImanager::GetUIpointer())->ApplyCommand("/tracking/verbose 0");
  }

  // Initialize user actions
  if(verbose > 0) {
    G4cout << "test31EventAction: Event # " 
           << nEvt << " started" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31EventAction::EndOfEventAction(const G4Event* evt)
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
    G4cout << "test31EventAction: Event # " 
           << nEvt << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  













