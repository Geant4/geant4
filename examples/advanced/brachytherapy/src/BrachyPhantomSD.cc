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
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ********************************
//    *                              *  
//    *    BrachyPhantomSD.cc       *
//    *                              *
//    ********************************
//
// $Id: BrachyPhantomSD.cc,v 1.5 2003/05/26 09:20:14 guatelli Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
#include "BrachyPhantomSD.hh"
#include "BrachyPhantomHit.hh"

#include "BrachyDetectorConstruction.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"

//....

BrachyPhantomSD::BrachyPhantomSD(G4String name, G4int NumVoxelX, G4int NumVoxelZ)
	:G4VSensitiveDetector(name),numberOfVoxelsX(NumVoxelX),numberOfVoxelsZ(NumVoxelZ)
{
 G4String HCname;
 collectionName.insert(HCname = "PhantomHitsCollection");
 voxelID = new G4int[numberOfVoxelsX * numberOfVoxelsZ * numberOfVoxelsZ];
  phantomHitsCollection = NULL;
}

BrachyPhantomSD::~BrachyPhantomSD()
{
 delete[] voxelID;
}

void BrachyPhantomSD::Initialize(G4HCofThisEvent*)
{
  G4int numberOfVoxelsY = 300;
  phantomHitsCollection = new BrachyPhantomHitsCollection(SensitiveDetectorName,collectionName[0]);

  for(G4int k=0;k<numberOfVoxelsZ;k++)
    for(G4int i=0;i<numberOfVoxelsX;i++)
      for(G4int j=0;j<numberOfVoxelsY;j++)
	voxelID[i+k*numberOfVoxelsX+j*numberOfVoxelsY*numberOfVoxelsY] = -1;
}

G4bool BrachyPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "PhantomPhys")
    return false;

  G4double energyDeposit = aStep->GetTotalEnergyDeposit();
  if(energyDeposit == 0.)
    return false;

  G4VPhysicalVolume* physVol = ROhist->GetVolume();
  //G4VPhysicalVolume* mothVol = ROhist->GetVolume(1);
  
  // Read Voxel indexes: i is the x index, k is the z index
  G4int k = ROhist->GetReplicaNumber(1);
  G4int i = ROhist->GetReplicaNumber(2);
  G4int j= ROhist->GetReplicaNumber();
  if(voxelID[i+k*numberOfVoxelsX+j*numberOfVoxelsX*numberOfVoxelsX] == -1)
    {
      BrachyPhantomHit* PhantomHit = new BrachyPhantomHit(physVol->GetLogicalVolume(),i,j,k);
      
      G4RotationMatrix rotM;
      if(physVol->GetObjectRotation())
	rotM = *(physVol->GetObjectRotation());
      
      PhantomHit->SetEdep(energyDeposit);
      PhantomHit->SetPos(physVol->GetTranslation());
      PhantomHit->SetRot(rotM);

      G4int VoxelID = phantomHitsCollection->insert(PhantomHit);
      voxelID[i+k*numberOfVoxelsX+j*numberOfVoxelsX*numberOfVoxelsX] = VoxelID - 1;
    }
  else
    (*phantomHitsCollection)
         [voxelID[i+k*numberOfVoxelsX+j*numberOfVoxelsX*numberOfVoxelsX]]->AddEdep(energyDeposit);

  return true;
}

void BrachyPhantomSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = GetCollectionID(0); 
    }
  HCE->AddHitsCollection(HCID,phantomHitsCollection);
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



