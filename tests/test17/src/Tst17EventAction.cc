// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17EventAction.cc,v 1.2 1999-12-15 14:54:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst17EventAction.hh"

#include "Tst17RunAction.hh"

#include "Tst17CalorHit.hh"
#include "Tst17EventActionMessenger.hh"

#include "g4rw/tvordvec.h"

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

Tst17EventAction::Tst17EventAction(Tst17RunAction* Tst17RA)
:calorimeterCollID(-1),eventMessenger(NULL),
 verboselevel(0),drawFlag("all"),runaction (Tst17RA)
{
  eventMessenger = new Tst17EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst17EventAction::~Tst17EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17EventAction::BeginOfEventAction(const G4Event* evt)
{  if(calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 

  if(verboselevel>1)
    G4cout << "<<< Event  " << evt->GetEventID() << " started." << G4endl;
  nstep = 0. ;
  nstepCharged = 0. ;
  nstepNeutral = 0. ;
  Nch = 0. ;
  Nne = 0. ;
  NE=0.;
  NP=0.;
  Transmitted=0.;
  Reflected  =0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst17EventAction::EndOfEventAction(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Tst17CalorHitsCollection* CHC = NULL;
  if (HCE)
      CHC = (Tst17CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
   {
    int n_hit = CHC->entries();
   // if(verboselevel==2)
   // G4cout << "     " << n_hit
   //      << " hits are stored in Tst17CalorHitsCollection." << G4endl;

    G4double totEAbs=0, totLAbs=0;
    for (int i=0;i<n_hit;i++)
      { totEAbs += (*CHC)[i]->GetEdepAbs(); 
        totLAbs += (*CHC)[i]->GetTrakAbs();
      }
  if(verboselevel==2)
    G4cout
       << "   Absorber: total energy: " << G4std::setw(7) << 
                             G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << G4std::setw(7) <<
                             G4BestUnit(totLAbs,"Length")
       << G4endl;           

   // count event, add deposits to the sum ...
    runaction->CountEvent() ;
    runaction->AddTrackLength(totLAbs) ;
    runaction->AddnStepsCharged(nstepCharged) ;
    runaction->AddnStepsNeutral(nstepNeutral) ;
    if(verboselevel==2)
      G4cout << " Ncharged=" << Nch << "  ,   Nneutral=" << Nne << G4endl;
    runaction->CountParticles(Nch,Nne);
    runaction->AddEP(NE,NP);
    runaction->AddTrRef(Transmitted,Reflected) ;
    runaction->AddEdeps(totEAbs) ;
    runaction->FillEn(totEAbs) ;
  }

  if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  
  if(pVVisManager){

    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
    for(G4int i=0; i<n_trajectories; i++){ 
      G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") trj->DrawTrajectory(50);
      else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	trj->DrawTrajectory(50); 
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Tst17EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

void Tst17EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}
void Tst17EventAction::CountStepsCharged()
{
  nstepCharged += 1. ;
}
void Tst17EventAction::CountStepsNeutral()
{
  nstepNeutral += 1. ;
}
void Tst17EventAction::AddCharged() 
{
  Nch += 1.;
}
void Tst17EventAction::AddNeutral() 
{
  Nne += 1.;
}
void Tst17EventAction::AddE() 
{
  NE += 1.;
}
void Tst17EventAction::AddP() 
{
  NP += 1.;
}
void Tst17EventAction::SetTr()
{
  Transmitted = 1.;
}
void Tst17EventAction::SetRef()
{
  Reflected   = 1.;
}
  















