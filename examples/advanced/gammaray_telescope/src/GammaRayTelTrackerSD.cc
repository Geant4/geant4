// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelTrackerSD.cc,v 1.1 2001-03-05 13:58:23 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelTrackerSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelTrackerSD.hh"

#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"

#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerSD::GammaRayTelTrackerSD(G4String name,
                                   GammaRayTelDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  G4int NbOfTKRTiles  =  Detector->GetNbOfTKRTiles();
  NbOfTKRStrips  = Detector->GetNbOfTKRStrips();
  NbOfTKRLayers  = Detector->GetNbOfTKRLayers();  
  NbOfTKRStrips = NbOfTKRStrips*NbOfTKRTiles;  
  
  ThitXID = new G4int[NbOfTKRStrips][30];
  ThitYID = new G4int[NbOfTKRStrips][30];
  collectionName.insert("TrackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerSD::~GammaRayTelTrackerSD()
{
  delete [] ThitXID;
  delete [] ThitYID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  TrackerCollection = new GammaRayTelTrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);
 
 for (G4int i=0;i<NbOfTKRStrips;i++)
    for (G4int j=0;j<NbOfTKRLayers;j++) 
      {
	ThitXID[i][j] = -1;
	ThitYID[i][j] = -1;
      };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{ 
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ((edep/keV == 0.)) return false;      
  
  G4int StripTotal = Detector->GetNbOfTKRStrips();
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

  G4int StripNumber = 0;
  G4VPhysicalVolume* strip = 0;
  strip = ROhist->GetVolume();
  G4String StripName = strip->GetName();
  StripNumber= strip->GetCopyNo();  

  ROhist->MoveUpHistory();
  G4VPhysicalVolume* tile = ROhist->GetVolume(); 
  G4int TileNumber = tile->GetCopyNo();  
  G4String TileName = tile->GetName();   
  
  G4int NTile = (TileNumber%TileTotal);  
  G4int j=0;
      
  for (j=0;j<TileTotal;j++)
    { 
      if(NTile==j) StripNumber += StripTotal*NTile;
    }  
  
  // G4cout << " Plane Number = " << PlaneNumber << " " << PlaneName << G4endl;
  // G4cout << StripName << " " << StripNumber << G4endl;       
  
  ROhist->MoveUpHistory();
  
  G4VPhysicalVolume* ROPlane = ROhist->GetVolume(); 
  G4int ROPlaneNumber = ROPlane->GetCopyNo();
  G4String ROPlaneName = ROPlane->GetName();   
    
  if (PlaneName == "TKRDetectorX" )
    // The hit is on an X silicon plane
    {
      // This is a new hit
      if (ThitXID[StripNumber][PlaneNumber]==-1)
	{       
	  GammaRayTelTrackerHit* TrackerHit = new GammaRayTelTrackerHit;
	  TrackerHit->SetPlaneType(1);
	  TrackerHit->AddSil(edep);
	  TrackerHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  TrackerHit->SetNSilPlane(PlaneNumber);
	  TrackerHit->SetNStrip(StripNumber);
	  ThitXID[StripNumber][PlaneNumber] = 
	    TrackerCollection->insert(TrackerHit) -1;
	}
      else // This is not new
	{
	  (*TrackerCollection)[ThitXID[StripNumber][PlaneNumber]]->AddSil(edep);
          // G4cout << "X" << PlaneNumber << " " << StripNumber << G4endl;
	}
    }
 
  if (PlaneName == "TKRDetectorY")
    // The hit is on an Y silicon plane    
    {   
      // This is a new hit
      if (ThitYID[StripNumber][PlaneNumber]==-1)
	{       
	  GammaRayTelTrackerHit* TrackerHit = new GammaRayTelTrackerHit;
	  TrackerHit->SetPlaneType(0);
	  TrackerHit->AddSil(edep);
	  TrackerHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  TrackerHit->SetNSilPlane(PlaneNumber);
	  TrackerHit->SetNStrip(StripNumber);
	  ThitYID[StripNumber][PlaneNumber] = 
	    TrackerCollection->insert(TrackerHit)-1;
	}
      else // This is not new
	{
	  (*TrackerCollection)[ThitYID[StripNumber][PlaneNumber]]->AddSil(edep);
          // G4cout << "Y" << PlaneNumber << " " << StripNumber << G4endl;
	}
    }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,TrackerCollection);


  for (G4int i=0;i<NbOfTKRStrips;i++) 
    for (G4int j=0;j<NbOfTKRLayers;j++)
      {
	ThitXID[i][j] = -1;
	ThitYID[i][j] = -1;
      };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....













