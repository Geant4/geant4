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
// $Id: HadrontherapyEventAction.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, G. Candiano, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// --------------------------------------------------------------
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4VVisManager.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyPhantomHit.hh"
#include "HadrontherapyPhantomSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyMatrix.hh"

HadrontherapyEventAction::HadrontherapyEventAction(HadrontherapyMatrix* matrixPointer) :
  drawFlag("all" )
{ 
  hitsCollectionID = -1;
  matrix = matrixPointer; 
}

HadrontherapyEventAction::~HadrontherapyEventAction()
{
}

void HadrontherapyEventAction::BeginOfEventAction(const G4Event* )
{
 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 if(hitsCollectionID == -1)
 	hitsCollectionID = pSDManager -> GetCollectionID("HadrontherapyPhantomHitsCollection");
}

void HadrontherapyEventAction::EndOfEventAction(const G4Event* evt)
{  
  if(hitsCollectionID < 0)
    return;

  G4HCofThisEvent* HCE = evt -> GetHCofThisEvent();
  HadrontherapyPhantomHitsCollection* CHC = NULL; 
 
  if(HCE)
    CHC = (HadrontherapyPhantomHitsCollection*)(HCE -> GetHC(hitsCollectionID));
  
  if(CHC)
    {
      if(matrix)
	{
	  // Fill the matrix with the information: voxel and associated energy deposit 
          // in the phantom at the end of the event

	  G4int HitCount = CHC -> entries();
	  for (G4int h=0; h<HitCount; h++)
	    {
	      G4int i = ((*CHC)[h]) -> GetXID();
	      G4int j = ((*CHC)[h]) -> GetYID();
	      G4int k = ((*CHC)[h]) -> GetZID();
              G4double energyDeposit = ((*CHC)[h]) -> GetEdep();
              matrix -> Fill(i, j, k, energyDeposit);              
	    }
	}
    }

  // Extract the trajectories and draw them in the visualisation

  if (G4VVisManager::GetConcreteInstance())
    {
      G4TrajectoryContainer * trajectoryContainer = evt -> GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer -> entries();
      
      for (G4int i = 0; i < n_trajectories; i++) 
        {
          G4Trajectory* trj = (G4Trajectory*)
	    ((*(evt -> GetTrajectoryContainer()))[i]);
	  if(drawFlag == "all") trj -> DrawTrajectory(50);
	  else if((drawFlag == "charged")&&(trj -> GetCharge() != 0.))
	    trj -> DrawTrajectory(50);
	  else if ((drawFlag == "neutral")&&(trj -> GetCharge() == 0.))
	    trj -> DrawTrajectory(50);	     	     
	}
    }
}

