// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestCalorimeterSD -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestEventAction.hh"

#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::hTestEventAction(hTestRunAction* run):
  theRun(run),
  verbose(0),
  nEvt(0),
  drawFlag("all")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::~hTestEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::BeginOfEventAction(const G4Event* evt)
{  
  verbose = theRun->GetVerbose();
  theRun->InitializeTuples();
  nEvt++;

  if(verboselevel > 0) {
    G4cout << "hTestEventAction: Event #" 
           << evt->GetEventID() << " started; nEvt = " 
           << nEvt << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::EndOfEventAction(const G4Event* evt)
{

  theRun->FillTuples();

  /*
    G4int i,j;
    G4double EE[60];
    G4String EEE[60]={
      "E0", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", 
      "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E19", 
      "E20", "E21", "E22", "E23", "E24", "E25", "E26", "E27", "E28", "E29", 
      "E30", "E31", "E32", "E33", "E34", "E35", "E36", "E37", "E38", "E39", 
      "E40", "E41", "E42", "E43", "E44", "E45", "E46", "E47", "E48", "E49", 
      "E50", "E51", "E52", "E53", "E54", "E55", "E56", "E57", "E58", "E59"};
    for(i=0; i<60; i++){
      EE[i]=0.0;
      for(j=0; j<i+1; j++) {
        EE[i] += EnergySlice[j];
      }  
      //      runaction->SaveToTuple(EEE[i],EE[i]);      
    }
    */

#ifdef G4VIS_USE  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    G4TrajectoryContainer* trjc = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trjc) n_trajectories = trjc->entries();  

    for(i=0; i<n_trajectories; i++) {
      G4Trajectory* t = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") t->DrawTrajectory(50);
      else if ((drawFlag == "charged")&&(t->GetCharge() != 0.))
                             t->DrawTrajectory(50); 
    }
  }  
#endif

  if(verboselevel > 0) {
    G4cout << "hTestEventAction: Event #" 
           << evt->GetEventID() << " ended" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  













