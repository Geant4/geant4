//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "fluoTestSampleSD.hh"
#include "fluoTestSensorHit.hh"
#include "fluoTestDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestSampleSD::fluoTestSampleSD(G4String name,
                                   fluoTestDetectorConstruction* det)
:G4VSensitiveDetector(name),Sample(det)
{
  collectionName.insert("SamCollection");
  HitSID = new G4int[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestSampleSD::~fluoTestSampleSD()
{
  delete [] HitSID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSampleSD::Initialize(G4HCofThisEvent*HCE)
{
  SamCollection = new fluoTestSensorHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 

 {HitSID[0] = -1;};
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool fluoTestSampleSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
      stepl = aStep->GetStepLength();
      
  if ((edep==0.)&&(stepl==0.)) return false; 

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  //theTouchable->MoveUpHistory();
 
  if (HitSID[0]==-1)
    { 
      fluoTestSensorHit* samHit = new fluoTestSensorHit();
      if (physVol == Sample->GetSample())
      samHit->AddSam(edep,stepl);
      HitSID[0] = SamCollection->insert(samHit) - 1;
      if (verboseLevel>0)
        G4cout << " New Sample Hit  " << G4endl;
    }
 
else
    { 
      if (physVol == Sample->GetSample())
	{(*SamCollection)[HitSID[0]]->AddSam(edep,stepl);}
      if (verboseLevel>0)
        {G4cout << " Energy added to the sample "  << G4endl;} 
    }
 
 return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSampleSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,SamCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSampleSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSampleSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestSampleSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

