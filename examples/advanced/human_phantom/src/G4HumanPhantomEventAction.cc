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
//
 
#include "G4HumanPhantomEventAction.hh"
#include "G4HumanPhantomAnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"


G4HumanPhantomEventAction::G4HumanPhantomEventAction():
  eventNumber(0)
{ }
 
G4HumanPhantomEventAction::~G4HumanPhantomEventAction()
{ }

void G4HumanPhantomEventAction::BeginOfEventAction(const G4Event* evt)
{
 eventNumber = evt -> GetEventID();
 path = 0.;
}
 
void G4HumanPhantomEventAction::EndOfEventAction(const G4Event* evt)
{  
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  if (G4VVisManager::GetConcreteInstance())
    {
     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                            ((*(evt->GetTrajectoryContainer()))[i]);
          trj->DrawTrajectory(50);
        }
    }

#ifdef G4ANALYSIS_USE
  G4HumanPhantomAnalysisManager* analysis = G4HumanPhantomAnalysisManager::getInstance();
  if (path != 0.)analysis -> particlePath(path);
#endif
}

void G4HumanPhantomEventAction::SetPath(G4double particlePath)
{
  path = path + particlePath;
}
