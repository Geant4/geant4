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
//
//
#include "MedLinacPhantomSD.hh"
#include "MedLinacPhantomHit.hh"

#include "MedLinacDetectorConstruction.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"

//....

MedLinacPhantomSD::MedLinacPhantomSD(G4String name, G4int NumVoxelX, G4int NumVoxelY, G4int NumVoxelZ)
	:G4VSensitiveDetector(name),NumberOfVoxelsX(NumVoxelX),NumberOfVoxelsY(NumVoxelY),NumberOfVoxelsZ(NumVoxelZ)
{
 G4String HCname;
 collectionName.insert(HCname = "PhantomHitsCollection");
 voxelID = new G4int[NumberOfVoxelsX * NumberOfVoxelsY * NumberOfVoxelsZ];
  phantomHitsCollection = NULL;
}

MedLinacPhantomSD::~MedLinacPhantomSD()
{
  delete[] voxelID;
}

void MedLinacPhantomSD::Initialize(G4HCofThisEvent*)
{

  phantomHitsCollection = new MedLinacPhantomHitsCollection(SensitiveDetectorName,collectionName[0]);

  for(G4int k=0;k<NumberOfVoxelsZ;k++)
    for(G4int i=0;i<NumberOfVoxelsX;i++)
      for(G4int j=0;j<NumberOfVoxelsY;j++)
	voxelID[i+k*NumberOfVoxelsX+j*NumberOfVoxelsY*NumberOfVoxelsY] = -1;
}

G4bool MedLinacPhantomSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  if(!ROhist)
    return false;

  if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "Phantom_phys")
    return false;

  G4double energyDeposit = aStep->GetTotalEnergyDeposit();
  G4cout << "energy deposit in PhantomSD e' "<< energyDeposit << G4endl;
  if(energyDeposit == 0.)
    return false;

  G4VPhysicalVolume* physVol = ROhist->GetVolume();
  //G4VPhysicalVolume* mothVol = ROhist->GetVolume(1);
  
  // Read Voxel indexes: i is the x index, k is the z index
  G4int k = ROhist->GetReplicaNumber(1);
  G4int i = ROhist->GetReplicaNumber(2);
  G4int j = ROhist->GetReplicaNumber();
  if(voxelID[i+k*NumberOfVoxelsX+j*NumberOfVoxelsX*NumberOfVoxelsX] == -1)
    {
      MedLinacPhantomHit* PhantomHit = new MedLinacPhantomHit(physVol->GetLogicalVolume(),i,j,k);
      
      G4RotationMatrix rotM;
      if(physVol->GetObjectRotation())
	rotM = *(physVol->GetObjectRotation());
      
      PhantomHit->SetEdep(energyDeposit);
  G4cout << "SetEdep in  PhantomSD e----------' "<< energyDeposit << G4endl;
      PhantomHit->SetPos(physVol->GetTranslation());
      PhantomHit->SetRot(rotM);

      G4int VoxelID = phantomHitsCollection->insert(PhantomHit);
      voxelID[i+k*NumberOfVoxelsX+j*NumberOfVoxelsX*NumberOfVoxelsX] = VoxelID - 1;
    }
  else
    (*phantomHitsCollection)
         [voxelID[i+k*NumberOfVoxelsX+j*NumberOfVoxelsX*NumberOfVoxelsX]]->AddEdep(energyDeposit);

  return true;
}

void MedLinacPhantomSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = GetCollectionID(0); 
    }
  HCE->AddHitsCollection(HCID,phantomHitsCollection);
}

void MedLinacPhantomSD::clear()
{
} 

void MedLinacPhantomSD::DrawAll()
{
}

void MedLinacPhantomSD::PrintAll()
{
}



