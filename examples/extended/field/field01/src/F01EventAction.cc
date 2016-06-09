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
// $Id: F01EventAction.cc,v 1.6 2006/06/29 17:16:29 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F01EventAction.hh"

#include "F01RunAction.hh"

#include "F01CalorHit.hh"
#include "F01EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01EventAction::F01EventAction(F01RunAction* F01RA)
 : calorimeterCollID(-1), eventMessenger(0),
   runaction(F01RA), verboselevel(0),
   drawFlag("all"), printModulo(10000)
{
  eventMessenger = new F01EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01EventAction::~F01EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01EventAction::BeginOfEventAction(const G4Event* evt)
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01EventAction::EndOfEventAction(const G4Event* evt)
{
  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for(G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory(50);
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(50); 
      }
  }  

  if(verboselevel>0)
    G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
  
  
  //save rndm status
  if (runaction->GetRndmFreq() == 2)
    { 
     CLHEP::HepRandom::saveEngineStatus("endOfEvent.rndm");   
     G4int evtNb = evt->GetEventID();
     if (evtNb%printModulo == 0)
       { 
         G4cout << "\n---> End of Event: " << evtNb << G4endl;
         CLHEP::HepRandom::showEngineStatus();
       }
    }     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int F01EventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01EventAction::setEventVerbose(G4int level)
{
  verboselevel = level ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
