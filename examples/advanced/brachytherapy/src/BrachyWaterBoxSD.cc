//    ********************************
//    *                              *  
//    *    BrachyWaterBoxSD.cc       *
//    *                              *
//    ********************************

#include "BrachyWaterBoxSD.hh"
#include "BrachyWaterBoxHit.hh"
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

BrachyWaterBoxSD::BrachyWaterBoxSD(G4String name, G4int NumVoxelX, G4int NumVoxelZ)
	:G4VSensitiveDetector(name),m_NumVoxelX(NumVoxelX),m_NumVoxelZ(NumVoxelZ)
{
 G4String HCname;
 collectionName.insert(HCname="WaterBoxHitsCollection");
 m_pVoxelID = new G4int[NumVoxelX*NumVoxelZ];
 m_pWaterBoxHitsCollection = NULL;
}

//....

BrachyWaterBoxSD::~BrachyWaterBoxSD()
{
 delete[] m_pVoxelID;
}

//....

void BrachyWaterBoxSD::Initialize(G4HCofThisEvent*HCE)
{
 m_pWaterBoxHitsCollection = new BrachyWaterBoxHitsCollection(SensitiveDetectorName,collectionName[0]);

 for(G4int k=0;k<m_NumVoxelZ;k++)
 	for(G4int i=0;i<m_NumVoxelX;i++)
 		m_pVoxelID[i+k*m_NumVoxelX] = -1;
}

//....

G4bool BrachyWaterBoxSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
 if(!ROhist)
	return false;

 if(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "WaterBoxPhys")
	return false;

 G4double edep = aStep->GetTotalEnergyDeposit();
 if(edep==0.)
	return false;

 G4VPhysicalVolume* physVol = ROhist->GetVolume();
 G4VPhysicalVolume* mothVol = ROhist->GetVolume(1);
  
 // Read Voxel indexes: i is the x index, k is the z index
 G4int k = ROhist->GetReplicaNumber();
 G4int i = ROhist->GetReplicaNumber(1);

 if(m_pVoxelID[i+k*m_NumVoxelX]==-1)
 	{
	BrachyWaterBoxHit* WaterBoxHit = new BrachyWaterBoxHit(physVol->GetLogicalVolume(),i,k);
      
	G4RotationMatrix rotM;
        if(physVol->GetObjectRotation())
		rotM = *(physVol->GetObjectRotation());
      
	WaterBoxHit->SetEdep(edep);
        WaterBoxHit->SetPos(physVol->GetTranslation());
        WaterBoxHit->SetRot(rotM);

        G4int VoxelID = m_pWaterBoxHitsCollection->insert(WaterBoxHit);
        m_pVoxelID[i+k*m_NumVoxelX] = VoxelID - 1;
	}
 else
	(*m_pWaterBoxHitsCollection)[m_pVoxelID[i+k*m_NumVoxelX]]->AddEdep(edep);

 return true;
}

//....

void BrachyWaterBoxSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  	{ 
	HCID = GetCollectionID(0); 
	}
  HCE->AddHitsCollection(HCID,m_pWaterBoxHitsCollection);
}

//....

void BrachyWaterBoxSD::clear()
{
} 

//....

void BrachyWaterBoxSD::DrawAll()
{
} 

//....

void BrachyWaterBoxSD::PrintAll()
{
}
