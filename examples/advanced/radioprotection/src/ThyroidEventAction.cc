//    *******************************
//    *                             *
//    *    ThyroidEventAction.cc    *
//    *                             *
//    *******************************

#include "ThyroidEventAction.hh"
#include "ThyroidHit.hh"
#include "ThyroidSD.hh"

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

ThyroidEventAction::ThyroidEventAction() 
       
{
 m_HitsCollectionID = -1;
 //m_pVoxel = pVoxel;
 // SDname=SDName;
 //pDetector=new ThyroidDetectorConstruction(SDname);

 /* m_NumVoxelX=pDetector-> GetNumVoxelX();
  m_NumVoxelY=pDetector-> GetNumVoxelY();
  m_NumVoxelZ=pDetector->GetNumVoxelZ();
  VoxelWidth_X= pDetector -> VoxelWidth_X();
  VoxelWidth_Y= pDetector -> VoxelWidth_Y();
  VoxelWidth_Z=pDetector->VoxelWidth_Z();*/
 
}

//....

ThyroidEventAction::~ThyroidEventAction()
{
  // delete pDetector;
}

//....

void ThyroidEventAction::BeginOfEventAction(const G4Event* )
{
 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 if(m_HitsCollectionID == -1)
 	m_HitsCollectionID =
pSDManager->GetCollectionID("RightThyroidHitsCollection"); }


//....

void ThyroidEventAction::EndOfEventAction(const G4Event* evt)
{
G4int event_id = evt->GetEventID();
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
 G4cout << "    " << n_trajectories 
    << " trajectories stored in this event." << G4endl;

  if (event_id < 100 || event_id%100 == 0)
  {
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  G4cout << "    " << n_trajectories 
     << " trajectories stored in this event." << G4endl;
   }

 if(m_HitsCollectionID < 0)
	return;

 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 ThyroidHitsCollection* CHC = NULL; 
 
 if(HCE)

	CHC = (ThyroidHitsCollection*)(HCE->GetHC(m_HitsCollectionID));

 if(CHC)
	{
	  //if(m_pVoxel)
		{
		// Fill voxel matrix with energy deposit data
		G4int HitCount = CHC->entries();
		for (G4int h=0; h<HitCount; h++)
 
          G4int NumVoxelX=200;
          G4int NumVoxelY=600;
          G4int NumVoxelZ=1200;

	  G4double VoxelWidth_X=0.05*mm;
          G4double VoxelWidth_Y=0.05*mm;
          G4double VoxelWidth_Z=0.05*mm;
          G4int i=((*CHC)[h])->GetXID();
          G4int  w=((*CHC)[h])->GetYID();
          G4int  k=((*CHC)[h])->GetZID();

          G4double EDep=(*CHC)[h]->GetEdep();
          G4double x=(- NumVoxelX+1+2*i)*VoxelWidth_X/2;
          G4double y=(- NumVoxelY+1+2*w)*VoxelWidth_Y/2; 
          G4double z=(- NumVoxelZ+1+2*k)*VoxelWidth_Z/2;

		  //m_pVoxel[((*CHC)[h])->GetXID() + ((*CHC)[h])->GetYID() + ((*CHC)[h])->GetZID()] += (*CHC)[h]->GetEdep();
 if(EDep!=0)
G4cout<<i<<" "<<w<<""<<k<<" "<<x<<" "<<y<<" "<<z<<" "<<  EDep<<"Energia dello step" << G4endl;}
		
	}	
}










