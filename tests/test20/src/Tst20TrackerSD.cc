// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20TrackerSD.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20TrackerSD  ------
//
// ************************************************************

#include "Tst20TrackerSD.hh"

#include "Tst20TrackerHit.hh"
#include "Tst20DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"

#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20TrackerSD::Tst20TrackerSD(G4String name,
                                   Tst20DetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  G4int NbOfTKRTiles  =  Detector->GetNbOfTKRTiles();
  NbOfTKRPixels  = Detector->GetNbOfTKRPixels();
  NbOfTKRLayers  = Detector->GetNbOfTKRLayers();  
  NbOfTKRPixels = NbOfTKRPixels*NbOfTKRTiles;  
  HitID = new G4int[NbOfTKRPixels][30];

  collectionName.insert("TrackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20TrackerSD::~Tst20TrackerSD()
{
  delete [] HitID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  TrackerCollection = new Tst20TrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  for (G4int i=0;i<NbOfTKRPixels;i++)
    for (G4int j=0;j<NbOfTKRLayers;j++) 
      {
	HitID[i][j] = -1;
      };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool Tst20TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{ 
  
  G4double edep = aStep->GetTotalEnergyDeposit();

  if ((edep/keV == 0.)) return false;      
  
  G4int PixelsTotal = Detector->GetNbOfTKRPixels();
  G4int TileTotal  = Detector->GetNbOfTKRTiles();  
 
  // This TouchableHistory is used to obtain the physical volume
  // of the hit

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
 
  G4VPhysicalVolume* phys_tile = theTouchable->GetVolume();  
  G4VPhysicalVolume* plane = phys_tile->GetMother();  
  
  G4int PlaneNumber = 0;
  PlaneNumber=plane->GetCopyNo();
  G4String PlaneName = plane->GetName();
  
  // The RO History is used to obtain the real strip
  // of the hit

  G4int PixelNumber = 0;
  G4VPhysicalVolume* pixel = 0;
  pixel = ROhist->GetVolume();
  G4String PixelName = pixel->GetName();
  PixelNumber= pixel->GetCopyNo();  
  
  ROhist->MoveUpHistory();
  G4VPhysicalVolume* tile = ROhist->GetVolume(); 
  G4int TileNumber = tile->GetCopyNo();  
  G4String TileName = tile->GetName();   
  
  G4int NTile = (TileNumber%TileTotal);       // Tile Row
  G4int NPixel = (PixelNumber%PixelsTotal);   // Pixel Row
  
  //  PixelNumber += ((TileNumber+NPixel+NTile*PixelsTotal)*PixelsTotal);
 

  // G4cout << " Plane Number = " << PlaneNumber << " " << PlaneName << G4endl;
  // G4cout << StripName << " " << StripNumber << G4endl;       
  
  ROhist->MoveUpHistory();
  
  G4VPhysicalVolume* ROPlane = ROhist->GetVolume(); 
  G4int ROPlaneNumber = ROPlane->GetCopyNo();
  G4String ROPlaneName = ROPlane->GetName();   
    
  if (PlaneName == "TKRDetector" )

    // The hit is on an silicon plane
    {
      // This is a new hit
      if (HitID[PixelNumber][PlaneNumber]==-1)
	{       
	  Tst20TrackerHit* TrackerHit = new Tst20TrackerHit;
	  TrackerHit->AddSil(edep);
	  TrackerHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  TrackerHit->SetNSilPlane(PlaneNumber);
	  TrackerHit->SetNPixel(PixelNumber);
	  HitID[PixelNumber][PlaneNumber] = 
	    TrackerCollection->insert(TrackerHit) -1;
	}
      else // This is not new
	{
	  (*TrackerCollection)[HitID[PixelNumber][PlaneNumber]]->AddSil(edep);
	}
    }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20TrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,TrackerCollection);

  
  for (G4int i=0;i<NbOfTKRLayers;i++) 
    for (G4int j=0;j<NbOfTKRPixels;j++)
      {
	HitID[i][j] = -1;
      };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void Tst20TrackerSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20TrackerSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20TrackerSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....













