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
//    **********************************
//    *                                *
//    *    RemSimEventAction.cc        *
//    *                                *
//    **********************************
//
//
// $Id: RemSimEventAction.cc,v 1.7 2004/05/27 12:31:31 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Author : Susanna Guatelli, guatelli@ge.infn.it
// 
#include "RemSimEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandEngine.h"
RemSimEventAction::RemSimEventAction():evtNo(-1)
{}
 
RemSimEventAction::~RemSimEventAction()
{}
 
void RemSimEventAction::BeginOfEventAction(const G4Event* evt)
{ 
  evtNo = evt -> GetEventID(); 
  G4int printModul = 500;
  if (evtNo%printModul == 0) 
   G4cout << "\n---> Begin Of Event: " << evtNo << G4endl;
}

void RemSimEventAction::EndOfEventAction(const G4Event* evt)
{
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories =0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  if (G4VVisManager::GetConcreteInstance())
    {
      for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	    ((*(evt->GetTrajectoryContainer()))[i]);
	trj->DrawTrajectory(50);
        }
    }
}
