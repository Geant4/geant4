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
//    *******************************
//    *                             *
//    *    BrachyEventAction.cc  	*
//    *                             *
//    *******************************

#include "BrachyEventAction.hh"
#include "BrachyWaterBoxHit.hh"
#include "BrachyWaterBoxSD.hh"
#include "BrachyDetectorConstruction.hh"
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
//#include"BrachyAnalysisManager.hh"
//....

BrachyEventAction::BrachyEventAction(G4String &SDName) :
 drawFlag("all" ),printModulo(1)      
{
 m_HitsCollectionID = -1;

 SDname=SDName;
 pDetector=new BrachyDetectorConstruction(SDname);
 
   m_NumVoxelX=pDetector-> GetNumVoxelX();
  m_NumVoxelZ=pDetector->GetNumVoxelZ();
  if (!m_pVoxel){
                m_pVoxel=new G4float[ m_NumVoxelX * m_NumVoxelZ];
                for (G4int i=0;i<m_NumVoxelX * m_NumVoxelZ;i++)
		  { m_pVoxel[i]=0;}
  }
}

//....

BrachyEventAction::~BrachyEventAction()
{
 delete pDetector;
 if(m_pVoxel)delete[]m_pVoxel;
}


G4int BrachyEventAction::GetEventno()
{
  G4int evno = fpEventManager->GetConstCurrentEvent()->GetEventID() ;
  return evno ;
}


G4double BrachyEventAction::GetEnergy(G4int k) const
{ return m_pVoxel[k];}



//....

void BrachyEventAction::BeginOfEventAction(const G4Event* aEvent)
{

 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 if(m_HitsCollectionID == -1)
 	m_HitsCollectionID = pSDManager->GetCollectionID("WaterBoxHitsCollection");
}

//....

void BrachyEventAction::EndOfEventAction(const G4Event* evt)
{
 if(m_HitsCollectionID < 0)
	return;

 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 BrachyWaterBoxHitsCollection* CHC = NULL; 
 
 if(HCE)
	CHC = (BrachyWaterBoxHitsCollection*)(HCE->GetHC(m_HitsCollectionID));

 if(CHC)
	{
	if(m_pVoxel)
		{
		// Fill voxel matrix with energy deposit data
		G4int HitCount = CHC->entries();
		for (G4int h=0; h<HitCount; h++)
                {
		j=((*CHC)[h])->GetZID() + ((*CHC)[h])->GetXID()*m_NumVoxelX;
                m_pVoxel[j]+=(*CHC)[h]->GetEdep();
		}
		}
     }



  // extract the trajectories and draw them

  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)((*(evt->
						 GetTrajectoryContainer()))[i]);

 
 if (drawFlag == "all") trj->DrawTrajectory(50);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(50);
          else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
                                  trj->DrawTrajectory(50);	     	     
	}
    }
 

}









