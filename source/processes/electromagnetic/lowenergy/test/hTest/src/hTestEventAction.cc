// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestEventAction.hh"

#include "hTestRunAction.hh"

#include "hTestCalorHit.hh"
#include "hTestEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::hTestEventAction(hTestRunAction* hTestRA)
:calorimeterCollID(-1),eventMessenger(NULL),
 verboselevel(0),runaction(hTestRA),drawFlag("all")
{
  eventMessenger = new hTestEventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestEventAction::~hTestEventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::BeginOfEventAction(const G4Event* evt)
{  

  EnergyDeposition = 0.0;

  if(verboselevel>1)
    G4cout << " Event  " << evt->GetEventID() << " started!!!" << G4endl;
  nstep = 0. ;
  nstepCharged = 0. ;
  nstepNeutral = 0. ;
  Nch = 0. ;
  Nne = 0. ;
  NE=0.;
  NP=0.;
  Transmitted=0.;
  Reflected  =0.;
  for(G4int i=0; i<60; i++){EnergySlice[i] = 0.0;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::EndOfEventAction(const G4Event* evt)
{
   // count event, add deposits to the sum ...
    runaction->CountEvent() ;
    runaction->AddnStepsNeutral(nstepNeutral) ;
    if(verboselevel==2)
      G4cout << " Ncharged=" << Nch << "  ,   Nneutral=" << Nne << G4endl;
    runaction->CountParticles(Nch,Nne);
    runaction->AddEP(NE,NP);
    runaction->AddTrRef(Transmitted,Reflected) ;

    //    runaction->FillEn(EnergyDeposition,EnergyDeposition);
    nstep=nstepCharged+nstepNeutral ;
    runaction->FillNbOfSteps(nstep);
    runaction->SaveToTuple("EDEP",EnergyDeposition);      
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

    runaction->SaveEvent();
    /* 
#ifdef G4VIS_USE  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for(i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory(50);
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(50); 
      }
  }  
#endif
*/
  if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int hTestEventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::CountStepsCharged()
{
  nstepCharged += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::CountStepsNeutral()
{
  nstepNeutral += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::AddCharged() 
{
  Nch += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::AddNeutral() 
{
  Nne += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::AddE(G4double En, G4int copyNo) 
{
  EnergyDeposition += En;
  if(copyNo > -1 && copyNo < 60) EnergySlice[copyNo] += En;
  NE += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::AddP() 
{
  NP += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::SetTr()
{
  Transmitted = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestEventAction::SetRef()
{
  Reflected   = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  













