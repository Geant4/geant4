// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17CalorimeterSD.hh"

#include "Test17CalorHit.hh"
#include "Test17DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
  
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17CalorimeterSD::Test17CalorimeterSD(G4String name,
                                   Test17DetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("CalCollection");
  HitID = new G4int[500];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17CalorimeterSD::~Test17CalorimeterSD()
{
  delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17CalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  CalCollection = new Test17CalorHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for (G4int j=0;j<1; j++) {HitID[j] = -1;};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool Test17CalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();
      
  if ((edep==0.)&&(stepl==0.)) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4LogicalVolume* logVol = physVol->GetLogicalVolume();
  //theTouchable->MoveUpHistory();
  G4int Test17Number = 0 ;

  if (HitID[Test17Number]==-1)
    { 
      Test17CalorHit* calHit = new Test17CalorHit();
      if (logVol == Detector->GetAbsorber()) calHit->AddAbs(edep,stepl);
      HitID[Test17Number] = CalCollection->insert(calHit) - 1;
      if (verboseLevel>0)
        G4cout << " New Calorimeter Hit on Test17: " << Test17Number << G4endl;
    }
  else
    { 
      if (logVol == Detector->GetAbsorber())
         (*CalCollection)[HitID[Test17Number]]->AddAbs(edep,stepl);
      if (verboseLevel>0)
        G4cout << " Energy added to Test17: " << Test17Number << G4endl; 
    }
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17CalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,CalCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17CalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Test17CalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







