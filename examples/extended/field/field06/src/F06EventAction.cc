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
/// \file field/field06/src/F06EventAction.cc
/// \brief Implementation of the F06EventAction class
//
//
//

#include "F06EventAction.hh"

#include "F06EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"

F06EventAction::F06EventAction()
 : verboselevel(0), printModulo(10), drawFlag("all")
{
  eventMessenger = new F06EventActionMessenger(this);
}

F06EventAction::~F06EventAction()
{
  delete eventMessenger;
}

void F06EventAction::BeginOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
     
 if(verboselevel>0)
    G4cout << "<<< Event  " << evtNb << " started." << G4endl;
}

void F06EventAction::EndOfEventAction(const G4Event* evt)
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();

   G4int n_trajectories = 0;

   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for(G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if(trj->GetParticleDefinition()->GetParticleName() != "neutron") continue;
        if (drawFlag == "all") trj->DrawTrajectory(50);
      }
  }  

  if (verboselevel>0)
     G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
  
}

G4int F06EventAction::GetEventNo()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID();
  return evno ;
}

void F06EventAction::SetEventVerbose(G4int level)
{
  verboselevel = level ;
}
