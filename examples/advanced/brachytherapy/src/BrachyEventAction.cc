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
//  S.Guatelli
//
//    *******************************
//    *                             *
//    *    BrachyEventAction.cc     *
//    *                             *
//    *******************************
//
// $Id: BrachyEventAction.cc,v 1.13 2002-11-27 11:11:22 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "BrachyPrimaryGeneratorActionI.hh"
#include "BrachyEventAction.hh"
#include "BrachyPhantomHit.hh"
#include "BrachyPhantomSD.hh"
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

#ifdef G4ANALYSIS_USE
#include"BrachyAnalysisManager.hh"
#endif

#include "BrachyPrimaryGeneratorActionIr.hh"
//....

BrachyEventAction::BrachyEventAction(G4String &SDName) :
  drawFlag("all" ),printModulo(1)      
{
 
  m_HitsCollectionID = -1;

 SDname=SDName;


 pDetector=new BrachyDetectorConstruction(SDname);

 
 m_NumVoxelX=pDetector-> GetNumVoxelX();
 m_NumVoxelZ=pDetector->GetNumVoxelZ();
 VoxelWidth_Z= 0.1*cm  ;
 VoxelWidth_X=0.1*cm;
 
}

//....

BrachyEventAction::~BrachyEventAction()
{
 delete pDetector;
 
}

//....

void BrachyEventAction::BeginOfEventAction(const G4Event* aEvent)
{

 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 if(m_HitsCollectionID == -1)
 	m_HitsCollectionID = pSDManager->GetCollectionID("PhantomHitsCollection");

}

//....

void BrachyEventAction::EndOfEventAction(const G4Event* evt)
{
 
if(m_HitsCollectionID < 0)
	return;

 G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
 BrachyPhantomHitsCollection* CHC = NULL; 
 
 if(HCE)
	CHC = (BrachyPhantomHitsCollection*)(HCE->GetHC(m_HitsCollectionID));

 if(CHC)
	{

	
	       G4int HitCount = CHC->entries();
		
   for (G4int h=0; h<HitCount; h++)
                {
	      
#ifdef G4ANALYSIS_USE			  
                  BrachyAnalysisManager* analysis = BrachyAnalysisManager::getInstance();	
#endif             
                     i=((*CHC)[h])->GetZID();
                     k=((*CHC)[h])->GetXID();
                     j=((*CHC)[h])->GetYID();  
                     

                      EnergyDep=((*CHC)[h]->GetEdep())/keV;
                      
                        x = (-m_NumVoxelZ+1+2*k)*VoxelWidth_X/2; 
                        z = (- m_NumVoxelZ+1+2*i)*VoxelWidth_Z/2;
                        y=(- m_NumVoxelZ+1+2*j)*VoxelWidth_Z/2;
#ifdef G4ANALYSIS_USE
			if(EnergyDep!=0)
                         
			  { 
                            
			    if (y<1.*mm){if (y> -1.*mm) 
                                  {analysis->hist(x,z,EnergyDep);}}
			  
			    
			  }
			 
			if(EnergyDep!=0)analysis->fill_Tuple(x,y,z,EnergyDep);
#endif 
		       
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










