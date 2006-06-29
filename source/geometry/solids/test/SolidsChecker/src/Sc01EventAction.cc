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
// $Id: Sc01EventAction.cc,v 1.3 2006-06-29 18:54:00 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              09  Sept. 1998 by S.Magni
// ------------------------------------------------------------

#include "Sc01EventAction.hh"
#include "Sc01EventActionMessenger.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include <iostream>

Sc01EventAction::Sc01EventAction()
  : fDrawFlag("all")
{
  eventMessenger = new Sc01EventActionMessenger(this);
}

Sc01EventAction::~Sc01EventAction()
{
  delete eventMessenger;
}

void Sc01EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 
 if ( evtNb % printModulo == 0 )
 { 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
 }
   //  G4cout << ">>> Start Event " << evt->GetEventID() << G4endl;
}

void Sc01EventAction::EndOfEventAction(const G4Event* evt)
{
  // G4cout << ">>> End Event " << evt->GetEventID() << G4endl;

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer) n_trajectories = trajectoryContainer->entries(); 

  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();

  if(pVisManager) 
  {
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;

    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

    for (G4int i=0; i<n_trajectories; i++)
    { 
      G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[i]);

      // if ( fDrawFlag == "all" )        
      trj->DrawTrajectory(50); 
                                  //     pVisManager->Draw(*trj,1000);
      /*
      else if ( (fDrawFlag == "charged")&&
                (trj->GetCharge() != 0.) )  pVisManager->Draw(*trj,1000);

      else if ((fDrawFlag == "neutral")&&
               (trj->GetCharge() == 0.))    pVisManager->Draw(*trj,1000);
      */
    }
  }		 
}





