#include "GammaRayTelPayloadSD.hh"

#include "GammaRayTelPayloadHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadSD::GammaRayTelPayloadSD(G4String name,
                                   GammaRayTelDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("PayloadCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadSD::~GammaRayTelPayloadSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::Initialize(G4HCofThisEvent*HCE)
{
  PayloadCollection = new GammaRayTelPayloadHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelPayloadSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  if ((edep==0.)) return false;      

  // This TouchableHistory is used to obtain the physical volume
  // of the hit
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

  // G4VPhysicalVolume* strip = ROhist->GetVolume(); 
  G4VPhysicalVolume* strip = theTouchable->GetVolume(); 
  
  G4int StripTotal = Detector->GetNbOfStrips();
  G4int TileTotal  = Detector->GetNbOfTKRTiles();

  G4int StripNumber = 0;
  G4int PlaneNumber = 0;
  StripNumber = theTouchable->GetReplicaNumber(1); // numero della strip
  G4String nome;
  nome = strip->GetName();
  G4cout << " Numero Strip = " << StripNumber <<G4endl;

  theTouchable->MoveUpHistory();
  // G4VPhysicalVolume* tile = ROhist->GetVolume(); // mattonella
  G4VPhysicalVolume* tile = theTouchable->GetVolume(); 
  G4int numero = tile->GetCopyNo(); // (da 0 a (NbOfTKRTiles*NbOfTKRTiles)-1)
  nome = tile->GetName();  
  G4int NumeroTile = (numero%TileTotal); // ritorna modulo  

  G4cout << " Numero Tile = " << NumeroTile <<G4endl;

  G4int j=0;
  for (j=0;j<TileTotal;j++)
    { 
      if(NumeroTile==j)StripNumber += StripTotal*NumeroTile;
    }  
  G4cout << " Numero Strip (2) = " << StripNumber <<G4endl;
  
  theTouchable->MoveUpHistory();
  // ROhist->MoveUpHistory();
  // G4VPhysicalVolume* plane = ROhist->GetVolume(); // piano
  G4VPhysicalVolume* plane = theTouchable->GetVolume(); 
  PlaneNumber=plane->GetCopyNo();
  nome = plane->GetName();
  G4cout << " Numero Piano = " << PlaneNumber <<G4endl;
     
  
  // The hit is on an X silicon plane
  if (nome == "TKRDetectorX" )
    { 
      GammaRayTelPayloadHit* PayloadHit = new GammaRayTelPayloadHit();
      PayloadHit->SetPlaneType(1);
      PayloadHit->AddSil(edep);
      PayloadHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
      PayloadHit->SetNStrip(StripNumber);
      PayloadHit->SetNSilPlane(PlaneNumber);
      PayloadCollection->insert(PayloadHit);
    }
  
  // The hit is on an Y silicon plane
  if (nome == "TKRDetectorY")
    { 
      GammaRayTelPayloadHit* PayloadHit = new GammaRayTelPayloadHit();
      PayloadHit->SetPlaneType(0);
      PayloadHit->AddSil(edep);
      PayloadHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
      PayloadHit->SetNStrip(StripNumber);
      PayloadHit->SetNSilPlane(PlaneNumber);
      PayloadCollection->insert(PayloadHit);
    }
  
    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,PayloadCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












