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
// Code developed by: S.Guatelli
//
//    ********************************
//    *                              *  
//    *    BrachyPhantomSD.cc       *
//    *                              *
//    ********************************
//
// $Id: BrachyPhantomSD.cc,v 1.3 2004/05/25 12:39:09 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
#include "BrachyPhantomSD.hh"
#include "BrachyPhantomHit.hh"
#include "BrachyAnalysisManager.hh"
#include "BrachyDetectorConstruction.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"

BrachyPhantomSD::BrachyPhantomSD(G4String name, G4int NumVoxelX, G4int NumVoxelZ)
	:G4VSensitiveDetector(name),numberOfVoxelsX(NumVoxelX),
         numberOfVoxelsZ(NumVoxelZ)
{
}

BrachyPhantomSD::~BrachyPhantomSD()
{
}

void BrachyPhantomSD::Initialize(G4HCofThisEvent*)
{
}

G4bool BrachyPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "PhantomPhys")
    return false;

  G4double energyDeposit = aStep -> GetTotalEnergyDeposit();
  if(energyDeposit == 0.)
    return false;
  
  if(energyDeposit != 0)                       
    {             
#ifdef G4ANALYSIS_USE
      // Read Voxel indexes: i is the x index, k is the z index
      G4int k = ROhist -> GetReplicaNumber(1);
      G4int i = ROhist -> GetReplicaNumber(2);
      G4int j = ROhist -> GetReplicaNumber();
  
      G4int numberOfVoxelZ = 300;
      G4double voxelWidthZ = 1. *mm;
          
      G4double x = (- numberOfVoxelZ+1+2*i)*voxelWidthZ/2; 
      G4double y = (- numberOfVoxelZ+1+2*j)*voxelWidthZ/2;
      G4double z = (- numberOfVoxelZ+1+2*k)*voxelWidthZ/2;
 
      BrachyAnalysisManager* analysis = 
	BrachyAnalysisManager::getInstance();   
      if (y < 0.8*mm && y > -0.8*mm)
	{ 
	  analysis -> FillHistogramWithEnergy(x,z,energyDeposit/MeV);
   
	  if (z < 0.8*mm && z > -0.8*mm)
	    analysis -> DoseDistribution(x,energyDeposit/MeV);
	}
#endif  
    }
  return true;
}

void BrachyPhantomSD::EndOfEvent(G4HCofThisEvent*)
{
}

void BrachyPhantomSD::clear()
{
} 

void BrachyPhantomSD::DrawAll()
{
}

void BrachyPhantomSD::PrintAll()
{
}



