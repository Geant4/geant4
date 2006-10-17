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
//    **************************************
//    *                                    *
//    *        CellEventAction.cc          *
//    *                                    *
//    **************************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------
 
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#ifdef G4ANALYSIS_USE
#include "CellAnalysisManager.hh"
#endif
#include "CellTrackerHit.hh"
#include "CellEventAction.hh"

CellEventAction::CellEventAction()
{ }

CellEventAction::~CellEventAction()
{ }

void CellEventAction::BeginOfEventAction(const G4Event*)
{ 
if (collisionID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      collisionID = SDman->GetCollectionID("TstCellCollection");
    }  

 totalEnergy = 0;
}
 
void CellEventAction::EndOfEventAction(const G4Event* evt)
{
  /*
 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  CellTrackerHitsCollection* hitCollection = 0;
  G4int hitNumber = 0;
    
  if (HCE) hitCollection = (CellTrackerHitsCollection*)
                                                     (HCE->GetHC(collisionID));
  if(hitCollection)
    {
      hitNumber = hitCollection->entries();
      for (G4int i=0;i<hitNumber;i++)
	{
	  totalEnergy += (*hitCollection)[i]->GetEdep();
	}
    }

  CellAnalysisManager* analysis = CellAnalysisManager::getInstance();
  analysis -> FillEnergyDeposit(totalEnergy/MeV);
  */
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  if (G4VVisManager::GetConcreteInstance())
    {
     for (G4int i=0; i<n_trajectories; i++) 
        { 
         G4Trajectory* trj = (G4Trajectory*)
	 ((*(evt->GetTrajectoryContainer()))[i]);
	 trj -> DrawTrajectory(50);
        }
    }
}

G4int CellEventAction::GetEventNo()
{
  G4int evno = fpEventManager -> GetConstCurrentEvent() -> GetEventID() ;
  return evno ;
}

