// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelPayloadSD.cc,v 1.5 2000-12-06 16:53:14 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelPayloadSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

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
  G4int NbOfTKRTiles  =  Detector->GetNbOfTKRTiles();
  NbOfTKRStrips  = Detector->GetNbOfTKRStrips();
  NbOfTKRLayers  = Detector->GetNbOfTKRLayers();  
  NbOfTKRStrips = NbOfTKRStrips*NbOfTKRTiles;  
  
  HitXID = new G4int[NbOfTKRStrips][30];
  HitYID = new G4int[NbOfTKRStrips][30];
  collectionName.insert("PayloadCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPayloadSD::~GammaRayTelPayloadSD()
{
  delete [] HitXID;
  delete [] HitYID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::Initialize(G4HCofThisEvent*HCE)
{
  PayloadCollection = new GammaRayTelPayloadHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  for (G4int i=0;i<NbOfTKRStrips;i++)
    for (G4int j=0;j<NbOfTKRLayers;j++) 
      {
	HitXID[i][j] = -1;
	HitYID[i][j] = -1;
      };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelPayloadSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
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
      if (HitXID[StripNumber][PlaneNumber]==-1)
	{       
	  GammaRayTelPayloadHit* PayloadHit = new GammaRayTelPayloadHit();
	  PayloadHit->SetPlaneType(1);
	  PayloadHit->AddSil(edep);
	  PayloadHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  PayloadHit->SetNSilPlane(PlaneNumber);
	  PayloadHit->SetNStrip(StripNumber);
	  HitXID[StripNumber][PlaneNumber] = 
	    PayloadCollection->insert(PayloadHit) -1;
	}
      else // This is not new
	{
	  (*PayloadCollection)[HitXID[StripNumber][PlaneNumber]]->AddSil(edep);
          // G4cout << "X" << PlaneNumber << " " << StripNumber << G4endl;
	}
    }
 
  if (PlaneName == "TKRDetectorY")
    // The hit is on an Y silicon plane    
    {   
      // This is a new hit
      if (HitYID[StripNumber][PlaneNumber]==-1)
	{       
	  GammaRayTelPayloadHit* PayloadHit = new GammaRayTelPayloadHit();
	  PayloadHit->SetPlaneType(0);
	  PayloadHit->AddSil(edep);
	  PayloadHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  PayloadHit->SetNSilPlane(PlaneNumber);
	  PayloadHit->SetNStrip(StripNumber);
	  HitYID[StripNumber][PlaneNumber] = 
	    PayloadCollection->insert(PayloadHit)-1;
	}
      else // This is not new
	{
	  (*PayloadCollection)[HitYID[StripNumber][PlaneNumber]]->AddSil(edep);
          // G4cout << "Y" << PlaneNumber << " " << StripNumber << G4endl;
	}
    }
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPayloadSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,PayloadCollection);


  for (G4int i=0;i<NbOfTKRLayers;i++) 
    for (G4int j=0;j<NbOfTKRStrips;j++)
      {
	HitXID[i][j] = -1;
	HitYID[i][j] = -1;
      };
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













