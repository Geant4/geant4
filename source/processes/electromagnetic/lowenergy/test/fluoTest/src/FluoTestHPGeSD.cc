//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestHPGeSD.hh"
#include "FluoTestSensorHit.hh"
#include "FluoTestDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestHPGeSD::FluoTestHPGeSD(G4String name,
                                   FluoTestDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[500];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestHPGeSD::~FluoTestHPGeSD()
{
  delete [] HitHPGeID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestHPGeSD::Initialize(G4HCofThisEvent*HCE)
 
//initializes HCE with the hits collection(s) created by this 
  //sensitive detector
{ 
  HPGeCollection = new FluoTestSensorHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
 for (G4int j=0;j<Detector->GetNbOfPixels();j++)
 {HitHPGeID [j]= -1;};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FluoTestHPGeSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{	
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ((edep==0.)) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  theTouchable->MoveUpHistory();     
  G4int PixelNumber = 0;
  if (Detector->GetNbOfPixels()>1) PixelNumber= physVol->GetCopyNo() ;
  if ( HitHPGeID[PixelNumber]==-1)
   { 
      FluoTestSensorHit* HPGeHit = new FluoTestSensorHit();
      HPGeHit->AddEnergy(edep);
     HitHPGeID[PixelNumber] = HPGeCollection->insert(HPGeHit) - 1;
      if (verboseLevel>0)
	G4cout << " New Hit on pixel: " << PixelNumber << G4endl;
 }
  else
    { 
 (*HPGeCollection)[HitHPGeID[PixelNumber]]->AddEnergy(edep);
 //G4double ED =(*HPGeCollection)[HitHPGeID[PixelNumber]]->GetEdepTot(); 
 if (verboseLevel>0)
	G4cout << " Energy added to Pixel: " << PixelNumber << G4endl; 
    }

 return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestHPGeSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,HPGeCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestHPGeSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestHPGeSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestHPGeSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




