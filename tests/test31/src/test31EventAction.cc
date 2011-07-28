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

void test31EventAction::BeginOfEventAction(const G4Event*)
{
  // New event
  test31Histo::GetPointer()->AddEvent();
  nEvt = test31Histo::GetPointer()->Event();

  //  if(nEvt == 9) G4cout << "Nan= " << std::sqrt(-9.0) << G4endl;

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
  if(verbose > 0 || nEvt/100*100 == nEvt) {
    G4cout << "test31EventAction: Event # "
           << nEvt << " started" << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31EventAction::EndOfEventAction(const G4Event*)
{
  if(verbose > 0) {
    G4cout << "test31EventAction: Event # " 
           << nEvt << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
