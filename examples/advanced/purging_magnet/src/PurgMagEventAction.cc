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
// Code developed by:
//  S.Larsson
//
//    ********************************
//    *                              *
//    *    PurgMagEventAction.cc     *
//    *                              *
//    ********************************
//
// $Id: PurgMagEventAction.cc,v 1.2 2004/06/18 09:17:57 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PurgMagEventAction.hh"

#include "PurgMagRunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"

#ifdef G4ANALYSIS_USE
#include"PurgMagAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagEventAction::PurgMagEventAction(PurgMagRunAction* run)
  :PurgMagRun(run),drawFlag("all"),printModulo(10000)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagEventAction::~PurgMagEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagEventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0) 
   G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagEventAction::EndOfEventAction(const G4Event* evt)
{  
  if (G4VVisManager::GetConcreteInstance())
    {
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
      
      for (G4int i=0; i<n_trajectories; i++) 
        { 
	  //G4cout<< "Iteration" << i <<G4endl;
	  G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
	  if (drawFlag == "all")
	    { 
	      trj->DrawTrajectory(50);
	    }
	  else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))	    
	    {
	      trj->DrawTrajectory(50);
	    }
	} 
      //save rndm status
      if (PurgMagRun->GetRndmFreq() == 2)
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
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

















