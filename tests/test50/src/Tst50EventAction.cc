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
// $Id: Tst50EventAction.cc,v 1.18 2003-07-03 13:43:10 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
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
#include "Tst50AnalysisManager.hh"
#endif
#include "Tst50TrackerHit.hh"
#include "Tst50EventAction.hh"

Tst50EventAction::Tst50EventAction()
{ }

Tst50EventAction::~Tst50EventAction()
{ }

void Tst50EventAction::BeginOfEventAction(const G4Event*)
{ 
if (collisionID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      collisionID = SDman->GetCollectionID("Tst50Collection");
    }  

 totalEnergy = 0;
}
 
void Tst50EventAction::EndOfEventAction(const G4Event* evt)
{
 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  
  Tst50TrackerHitsCollection* hitCollection = 0;
  G4int hitNumber = 0;
    
  if (HCE) hitCollection = (Tst50TrackerHitsCollection*)
                                                     (HCE->GetHC(collisionID));
  if(hitCollection)
    {
      hitNumber = hitCollection->entries();
      for (G4int i=0;i<hitNumber;i++)
	{
	  totalEnergy += (*hitCollection)[i]->GetEdep();
	}
    }

  Tst50AnalysisManager* analysis = Tst50AnalysisManager::getInstance();
  analysis -> FillEnergyDeposit(totalEnergy/MeV);

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

G4int Tst50EventAction::GetEventNo()
{
  G4int evno = fpEventManager -> GetConstCurrentEvent() -> GetEventID() ;
  return evno ;
}

