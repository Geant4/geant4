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

 
#include "MedLinacEventAction.hh"
#include "MedLinacPhantomHit.hh"
#include "MedLinacPhantomSD.hh"
#include "MedLinacDetectorConstruction.hh"
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

#ifdef G4ANALYSIS_USE
#include"MedLinacAnalysisManager.hh"
#endif


 // Retrieve information about the energy deposit in the phantom ...

MedLinacEventAction::MedLinacEventAction(G4String SDName) :
  drawFlag("all" )
{
  hitsCollectionID = -1;

  G4String sensitiveDetectorName=SDName;
  detector = MedLinacDetectorConstruction::GetInstance(sensitiveDetectorName);
  numberOfVoxelZ = detector->GetNumVoxelZ();
  voxelWidthZ = 0.1*cm; 
 }

 
MedLinacEventAction::~MedLinacEventAction()
{
  delete detector;
}


 
void MedLinacEventAction::BeginOfEventAction(const G4Event*)
{
  

  G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();
  if(hitsCollectionID == -1)
     hitsCollectionID = 
           sensitiveDetectorManager->GetCollectionID("PhantomHitsCollection");
  G4cout << " begin of event action hitsCollectionID "<< hitsCollectionID << G4endl;
}

 
void MedLinacEventAction::EndOfEventAction(const G4Event* evt)
{  
  if(hitsCollectionID < 0) return; 

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  MedLinacPhantomHitsCollection* CHC = NULL; 
  //G4cout << "prima dell'analisi HCE vale: "<< HCE  <<G4endl;
  if(HCE)
     CHC = (MedLinacPhantomHitsCollection*)(HCE->GetHC(hitsCollectionID));

  if(CHC)
    {
      G4int hitCount = CHC->entries();
      G4cout << " hitCount vale:"<< hitCount <<"----------------" <<G4endl;
      for (G4int h = 0; h < hitCount; h++)
	{

#ifdef G4ANALYSIS_USE	  
          //Store information about energy deposit
	  MedLinacAnalysisManager* analysis = 
                                      MedLinacAnalysisManager::getInstance();   
	  
          G4int i=((*CHC)[h])->GetZID();
	  G4int k=((*CHC)[h])->GetXID();
	  G4int j=((*CHC)[h])->GetYID();                       

	  G4double EnergyDep=((*CHC)[h]->GetEdep());
          G4cout << " EnergyDep in MedLinacEventAction e'-------------------- " << EnergyDep <<G4endl;       
	  G4double x = (-numberOfVoxelZ+1+2*k)*voxelWidthZ/2; 
	  G4double z = (- numberOfVoxelZ+1+2*i)*voxelWidthZ/2;
	  G4double y = (- numberOfVoxelZ+1+2*j)*voxelWidthZ/2;


	  //--------------------------------------------------------------
	     //2Dhistogram with the distribution of energy in the surface of the phantom 
	     //(x,y,energy)  (YThickness = 1. mm)
	  if(EnergyDep != 0)                       
	    {  
	      if (z<150.*mm){if (z> 149.*mm) 
		{analysis->FillHistogram3WithEnergy(x,y,EnergyDep/MeV);}
		}
	      }
	  //---------------------------------------------------------------
	  if(EnergyDep != 0)                       
	  {
	  analysis->FillHistogram1WithEnergy(x,z,EnergyDep/MeV);
	  }
	  //**** PDD in isocenter (Y and X Thickness = 1. mm) ***********			
 if(EnergyDep != 0)                       
	    { 
	      if (y<1.0*mm){if (y> -1.0*mm)
		{
		  if(x<1.0*mm){if (x> -1.0*mm)
		{analysis->FillHistogram4WithEnergy(z,EnergyDep/MeV);}
		  }}}
	    }
          //**************************************************************

          //***** flatness  along x **************************************

  if(EnergyDep != 0)                       
    {  
       if (z<150.*mm){if (z> 149.*mm) 
	 { if (y<1.0*mm){if (y> -1.0*mm)
       {analysis->FillHistogram5WithEnergy(x,EnergyDep/MeV);}

	                }
	 }
                     }
    }
          //**************************************************************
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

