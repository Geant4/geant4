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
// $Id: BrachyEventAction.cc,v 1.17 2003/05/27 08:37:54 guatelli Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
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

// Retrieve information about the energy deposit in the phantom ...

BrachyEventAction::BrachyEventAction(G4String &SDName) :
  drawFlag("all" )
{
  hitsCollectionID = -1;

  G4String sensitiveDetectorName = SDName;
  
  detector=new BrachyDetectorConstruction(sensitiveDetectorName);
  numberOfVoxelZ = detector->GetNumVoxelZ();
  voxelWidthZ = 0.1*cm;  
}

BrachyEventAction::~BrachyEventAction()
{
  delete detector;
}

void BrachyEventAction::BeginOfEventAction(const G4Event*)
{
  G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();
  if(hitsCollectionID == -1)
     hitsCollectionID = 
           sensitiveDetectorManager->GetCollectionID("PhantomHitsCollection");
}

void BrachyEventAction::EndOfEventAction(const G4Event* evt)
{  
  if(hitsCollectionID < 0) return; 

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  BrachyPhantomHitsCollection* CHC = NULL; 
 
  if(HCE)
     CHC = (BrachyPhantomHitsCollection*)(HCE->GetHC(hitsCollectionID));

  if(CHC)
    {
      G4int hitCount = CHC->entries();
      for (G4int h = 0; h < hitCount; h++)
	{
#ifdef G4ANALYSIS_USE	  
          //Store information about energy deposit in a 2DHistogram and in
	  // a ntuple ... 
	  BrachyAnalysisManager* analysis = 
                                      BrachyAnalysisManager::getInstance();   
	  
          G4int i=((*CHC)[h])->GetZID();
	  G4int k=((*CHC)[h])->GetXID();
	  G4int j=((*CHC)[h])->GetYID();                       

	  G4double EnergyDep=((*CHC)[h]->GetEdep());
                      
	  G4double x = (-numberOfVoxelZ+1+2*k)*voxelWidthZ/2; 
	  G4double z = (- numberOfVoxelZ+1+2*i)*voxelWidthZ/2;
	  G4double y = (- numberOfVoxelZ+1+2*j)*voxelWidthZ/2;

	  if(EnergyDep != 0)                       
	    { 
             if (y<1.*mm){if (y> -1.*mm) 
	      {analysis->FillHistogramWithEnergy(x,z,EnergyDep/MeV);}}
	     }
			 
	  if(EnergyDep != 0)analysis->FillNtupleWithEnergy(x,y,z,EnergyDep/MeV);
#endif 	       
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
