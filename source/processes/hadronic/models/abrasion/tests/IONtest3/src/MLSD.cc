////////////////////////////////////////////////////////////////////////////////
//
#include "MLSD.hh"

#include "MLGeometryConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
//#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLSD::MLSD (G4String name, MLGeometryConstruction* det)
  :G4VSensitiveDetector(name),geometry(det)
{
  collectionName.insert("MLCollection");
}
////////////////////////////////////////////////////////////////////////////////
//
MLSD::~MLSD ()
{
}
////////////////////////////////////////////////////////////////////////////////
//
void MLSD::Initialize (G4HCofThisEvent*HCE)
{
  static int HCID = -1;
  MLCollection = new MLHitsCollection(SensitiveDetectorName,collectionName[0]);

  if (HCID < 0) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,MLCollection);
}
////////////////////////////////////////////////////////////////////////////////
//
G4bool MLSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
  G4double edep = aStep->GetTotalEnergyDeposit() ;
  //                  aStep->GetTrack()->GetWeight();
  if (edep == 0.)  return false;
  const G4VTouchable* theTouchable = aStep->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume* physVol       = theTouchable->GetVolume();
  G4int LayerNumber                = physVol->GetCopyNo();
  G4double weight                  = aStep->GetTrack()->GetWeight();
  MLHit* aHit                      = new MLHit();
  aHit->SetEdep(LayerNumber,edep, weight);
  MLCollection->insert(aHit); 
  return true;
}
////////////////////////////////////////////////////////////////////////////////
