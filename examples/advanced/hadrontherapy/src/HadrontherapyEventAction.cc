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
// $Id: HadrontherapyEventAction.cc,v 2.0
// ------------------------------------------------------------
//     GEANT 4 - Hadrontherapy example
// ------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// -------------------------------------------------------------


#include "HadrontherapyEventAction.hh"
#include "HadrontherapyPhantomHit.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"

HadrontherapyEventAction::HadrontherapyEventAction(G4double *Matrix, G4int nbVoxelX, G4int nbVoxelY, G4int nbVoxelZ) :
  drawFlag("all" )
{ 
 m_HitsCollectionID = -1;
 
 numberX = nbVoxelX;
 numberY = nbVoxelY;
 numberZ = nbVoxelZ;
 matrix = Matrix; 
}

HadrontherapyEventAction::~HadrontherapyEventAction()
{
}

void HadrontherapyEventAction::BeginOfEventAction(const G4Event*evt )
{
  //G4int event_id = evt->GetEventID();
  //G4cout << event_id << ": evt" << G4endl;
 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 if(m_HitsCollectionID == -1)
 	m_HitsCollectionID = pSDManager->GetCollectionID("HadrontherapyPhantomHitsCollection");
}

void HadrontherapyEventAction::EndOfEventAction(const G4Event* evt)
{  
 if(m_HitsCollectionID < 0)
	return;

 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 HadrontherapyPhantomHitsCollection* CHC = NULL; 
 
 if(HCE)
	CHC = (HadrontherapyPhantomHitsCollection*)(HCE->GetHC(m_HitsCollectionID));

 if(CHC)
	{
	if(matrix)
		{
		// Fill voxel matrix with energy deposit data
		G4int HitCount = CHC->entries();
		for (G4int h=0; h<HitCount; h++)
		  {
		    G4int i = ((*CHC)[h])->GetXID();
                    G4int j = ((*CHC)[h])->GetYID();
                    G4int k =  ((*CHC)[h])->GetZID();
			matrix[(i * numberY + j)* numberZ + k]+= 
                        (*CHC)[h]->GetEdep();
			//	G4cout<<"Hit:"<< h << G4endl;
			//G4cout<< "Energy deposit in the event:" 
                        //      << matrix[(i * numberY + j)* numberZ + k]
			//     <<"in"<< i <<" " << j << " "<< "" << k<< G4endl;
		  }
		}
	}

  // extract the trajectories and draw them ...

  if (G4VVisManager::GetConcreteInstance())
    {
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
      
      for (G4int i=0; i<n_trajectories; i++) 
        {
          G4Trajectory* trj = (G4Trajectory*)
	    ((*(evt->GetTrajectoryContainer()))[i]);
	  if(drawFlag == "all") trj->DrawTrajectory(50);
	  else if((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	    trj->DrawTrajectory(50);
	  else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
	    trj->DrawTrajectory(50);	     	     
	}
    }

}

