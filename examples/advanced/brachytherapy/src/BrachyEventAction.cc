//    *******************************
//    *                             *
//    *    BrachyEventAction.cc  	*
//    *                             *
//    *******************************

#include "BrachyEventAction.hh"
#include "BrachyWaterBoxHit.hh"
#include "BrachyWaterBoxSD.hh"

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

//....

BrachyEventAction::BrachyEventAction(G4float *pVoxel,G4int NumVoxelX,G4int NumVoxelZ) :
	m_NumVoxelX(NumVoxelX),m_NumVoxelZ(NumVoxelZ)
{
 m_HitsCollectionID = -1;
 m_pVoxel = pVoxel;
}

//....

BrachyEventAction::~BrachyEventAction()
{
}

//....

void BrachyEventAction::BeginOfEventAction(const G4Event*)
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
			m_pVoxel[((*CHC)[h])->GetZID() + ((*CHC)[h])->GetXID()*m_NumVoxelX] += (*CHC)[h]->GetEdep();
		}
	}	
}
