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
//
// $Id: F02EventAction.cc,v 1.6 2006-06-29 18:27:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F02EventAction.hh"

#include "F02RunAction.hh"

#include "F02CalorHit.hh"
#include "F02EventActionMessenger.hh"


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
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02EventAction::F02EventAction(F02RunAction* F02RA)
:calorimeterCollID(-1),eventMessenger(NULL),
 verboselevel(0),runaction(F02RA),drawFlag("all"),printModulo(10000)
{
  eventMessenger = new F02EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F02EventAction::~F02EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
     
  if(verboselevel>1)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
    
  if (calorimeterCollID==-1)
    {
     G4SDManager * SDman = G4SDManager::GetSDMpointer();
     calorimeterCollID = SDman->GetCollectionID("CalCollection");
    } 

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

void F02EventAction::EndOfEventAction(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  F02CalorHitsCollection* CHC = NULL;
  
  if (HCE) CHC = (F02CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
  {
    int n_hit = CHC->entries();
   // if(verboselevel==2)
   // G4cout << "     " << n_hit
   //      << " hits are stored in F02CalorHitsCollection." << G4endl;

    G4double totEAbs=0, totLAbs=0;

    for (int i=0;i<n_hit;i++)
    { 
      totEAbs += (*CHC)[i]->GetEdepAbs(); 
      totLAbs += (*CHC)[i]->GetTrakAbs();
    }
    if(verboselevel==2)
    {
      G4cout << "   Absorber: total energy: " << std::setw(7) 
             << G4BestUnit(totEAbs,"Energy")
             << "       total track length: " << std::setw(7) 
             << G4BestUnit(totLAbs,"Length")  << G4endl;           
    }
    // count event, add deposits to the sum ...

    runaction->CountEvent() ;
    runaction->AddTrackLength(totLAbs) ;
    runaction->AddnStepsCharged(nstepCharged) ;
    runaction->AddnStepsNeutral(nstepNeutral) ;

    if(verboselevel==2)
    {
      G4cout << " Ncharged=" << Nch << "  ,   Nneutral=" << Nne << G4endl;
    }
    runaction->CountParticles(Nch,Nne);
    runaction->AddEP(NE,NP);
    runaction->AddTrRef(Transmitted,Reflected) ;
    runaction->AddEdeps(totEAbs) ;
    runaction->FillEn(totEAbs) ;

    nstep=nstepCharged+nstepNeutral ;
    runaction->FillNbOfSteps(nstep);
  }
  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;

    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries(); 
 
    for(G4int i=0; i<n_trajectories; i++) 
    { 
      G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);

      if (drawFlag == "all") 
      {
        trj->DrawTrajectory(50);
      }
      else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
      {
        trj->DrawTrajectory(50); 
      }
    }
  }  
  if(verboselevel>0)
  {  
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
  }
  
  //save rndm status
  if (runaction->GetRndmFreq() == 2)
  { 
    HepRandom::saveEngineStatus("endOfEvent.rndm");   
    G4int evtNb = evt->GetEventID();

    if (evtNb%printModulo == 0)
    { 
        G4cout << "\n---> End of Event: " << evtNb << G4endl;
        HepRandom::showEngineStatus();
    }
  }     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int F02EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::CountStepsCharged()
{
  nstepCharged += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::CountStepsNeutral()
{
  nstepNeutral += 1. ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::AddCharged() 
{
  Nch += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::AddNeutral() 
{
  Nne += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::AddE() 
{
  NE += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::AddP() 
{
  NP += 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::SetTr()
{
  Transmitted = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F02EventAction::SetRef()
{
  Reflected   = 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  











