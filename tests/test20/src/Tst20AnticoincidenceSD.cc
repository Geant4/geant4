// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20AnticoincidenceSD.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20AnticoincidenceSD  ------
//           by  G.Depaola & F.Longo (13 mar 2001)
//
// ************************************************************

#include "Tst20AnticoincidenceSD.hh"
#include "Tst20AnticoincidenceHit.hh"
#include "Tst20DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20AnticoincidenceSD::Tst20AnticoincidenceSD(G4String name,
                                   Tst20DetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  NbOfACDLateralTiles  =  Detector->GetNbOfACDLateralTiles();
  NbOfACDTopTiles  =  Detector->GetNbOfACDTopTiles(); 

  //G4cout <<  NbOfACDLateralTiles << " LAT " << G4endl;
  //G4cout <<  NbOfACDTopTiles << " TOP " << G4endl;
  
  HitLateralID = new G4int[NbOfACDLateralTiles];
  HitTopID = new G4int[NbOfACDTopTiles];
  collectionName.insert("AnticoincidenceCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20AnticoincidenceSD::~Tst20AnticoincidenceSD()
{
  delete [] HitLateralID;
  delete [] HitTopID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceSD::Initialize(G4HCofThisEvent*HCE)
{
  AnticoincidenceCollection = new Tst20AnticoincidenceHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  for (G4int i=0;i<NbOfACDLateralTiles;i++)
    {
      HitLateralID[i] = -1;
    };
  
  for (G4int j=0;j<NbOfACDTopTiles;j++) 
    {
      
      HitTopID[j] = -1;
    };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool Tst20AnticoincidenceSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{ 
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ((edep/keV == 0.)) return false;      
  
  // This TouchableHistory is used to obtain the physical volume
  // of the hit
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  
  G4VPhysicalVolume* acd_tile = theTouchable->GetVolume();  

  G4int ACDTileNumber=acd_tile->GetCopyNo();
  G4String ACDTileName = acd_tile->GetName();
  
  if (ACDTileName == "ACT" )
    // The hit is on an top ACD tile (ACDType 0)
    {
      // This is a new hit
      if (HitTopID[ACDTileNumber]==-1)
	{       
	  Tst20AnticoincidenceHit* 
	    AnticoincidenceHit = new Tst20AnticoincidenceHit;
	  AnticoincidenceHit->SetACDType(0);
	  AnticoincidenceHit->AddEnergy(edep);
	  AnticoincidenceHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  AnticoincidenceHit->SetACDTileNumber(ACDTileNumber);
	  HitTopID[ACDTileNumber] = 
	    AnticoincidenceCollection->insert(AnticoincidenceHit) -1;
	}
      else // This is not new
	{
	  (*AnticoincidenceCollection)
	    [HitTopID[ACDTileNumber]]->AddEnergy(edep);
	}
    }
  if (ACDTileName == "ACB" )
    // The hit is on an botton ACD tile (ACDType 0)
    {
      // This is a new hit
      if (HitTopID[ACDTileNumber]==-1)
	{       
	  Tst20AnticoincidenceHit* 
	    AnticoincidenceHit = new Tst20AnticoincidenceHit;
	  AnticoincidenceHit->SetACDType(0);
	  AnticoincidenceHit->AddEnergy(edep);
	  AnticoincidenceHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  AnticoincidenceHit->SetACDTileNumber(ACDTileNumber);
	  HitTopID[ACDTileNumber] = 
	    AnticoincidenceCollection->insert(AnticoincidenceHit) -1;
	}
      else // This is not new
	{
	  (*AnticoincidenceCollection)
	    [HitTopID[ACDTileNumber]]->AddEnergy(edep);
	}
    }
 
  if (ACDTileName == "ACL1")
    // The hit is on an lateral (left-right) ACD tile (ACDType 1)
    {   
      // This is a new hit
      if (HitLateralID[ACDTileNumber]==-1)
	{       
	  Tst20AnticoincidenceHit* 
	    AnticoincidenceHit = new Tst20AnticoincidenceHit;
	  AnticoincidenceHit->SetACDType(1);
	  AnticoincidenceHit->AddEnergy(edep);
	  AnticoincidenceHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  AnticoincidenceHit->SetACDTileNumber(ACDTileNumber);
	  HitLateralID[ACDTileNumber] = 
	    AnticoincidenceCollection->insert(AnticoincidenceHit) -1;
	}
      else // This is not new
	{
	  (*AnticoincidenceCollection)
	    [HitLateralID[ACDTileNumber]]->AddEnergy(edep);
	}
    }

   if (ACDTileName == "ACL2")
    // The hit is on an lateral (rear - front) ACD tile (ACDType 2)
    {   
      // This is a new hit
      if (HitLateralID[ACDTileNumber]==-1)
	{       
	  Tst20AnticoincidenceHit* 
	    AnticoincidenceHit = new Tst20AnticoincidenceHit;
	  AnticoincidenceHit->SetACDType(2);
	  AnticoincidenceHit->AddEnergy(edep);
	  AnticoincidenceHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  AnticoincidenceHit->SetACDTileNumber(ACDTileNumber);
	  HitLateralID[ACDTileNumber] = 
	    AnticoincidenceCollection->insert(AnticoincidenceHit) -1;
	}
      else // This is not new
	{
	  (*AnticoincidenceCollection)
	    [HitLateralID[ACDTileNumber]]->AddEnergy(edep);
	}
    }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,AnticoincidenceCollection);

  for (G4int i=0;i<NbOfACDLateralTiles;i++)
    {
      HitLateralID[i] = -1;
    };
  
  for (G4int j=0;j<NbOfACDTopTiles;j++) 
    {
      
      HitTopID[j] = -1;
    };

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20AnticoincidenceSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














