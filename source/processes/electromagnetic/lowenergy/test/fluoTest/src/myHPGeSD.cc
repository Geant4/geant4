//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "myHPGeSD.hh"
#include "mySensorHit.hh"
#include "myDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

myHPGeSD::myHPGeSD(G4String name,
                                   myDetectorConstruction* det)
:G4VSensitiveDetector(name),HPGe(det)
{
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[1];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

myHPGeSD::~myHPGeSD()
{
  delete [] HitHPGeID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myHPGeSD::Initialize(G4HCofThisEvent*HCE)
  
//initializes HCE with the hits collection(s) created by this 
  //sensitive detector
{
  HPGeCollection = new mySensorHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 

 {HitHPGeID[0] = -1;};
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool myHPGeSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
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
 
  if (HitHPGeID[0]==-1)
    { 
      mySensorHit* HPGeHit = new mySensorHit();
      if (physVol == HPGe->GetHPGe()) 
	//verify whether the physical volume is the detector's one
      HPGeHit->AddHPGe(edep,stepl);
      HitHPGeID[0] = HPGeCollection->insert(HPGeHit) - 1;
      if (verboseLevel>0)
        G4cout << " New HPGe Hit  " << G4endl;
    }
 
else
    { 
      if (physVol == HPGe->GetHPGe())
	{(*HPGeCollection)[HitHPGeID[0]]->AddHPGe(edep,stepl);}
      if (verboseLevel>0)
        {G4cout << " Energy added to the HPGe detector "  << G4endl;} 
    }
 
 return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myHPGeSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,HPGeCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myHPGeSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myHPGeSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myHPGeSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




