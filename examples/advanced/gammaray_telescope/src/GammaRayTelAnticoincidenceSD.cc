//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelAnticoincidenceSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
#include "G4RunManager.hh"
#include "GammaRayTelAnticoincidenceSD.hh"
#include "GammaRayTelAnticoincidenceHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceSD::GammaRayTelAnticoincidenceSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  Detector =
    (GammaRayTelDetectorConstruction*)(runManager->GetUserDetectorConstruction());
  
  NbOfACDLateralTiles  =  Detector->GetNbOfACDLateralTiles();
  NbOfACDTopTiles  =  Detector->GetNbOfACDTopTiles(); 

  G4cout <<  NbOfACDLateralTiles << " LAT " << G4endl;
  G4cout <<  NbOfACDTopTiles << " TOP " << G4endl;
  
  HitLateralID = new G4int[NbOfACDLateralTiles];
  HitTopID = new G4int[NbOfACDTopTiles];
  collectionName.insert("AnticoincidenceCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceSD::~GammaRayTelAnticoincidenceSD()
{
  delete [] HitLateralID;
  delete [] HitTopID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::Initialize(G4HCofThisEvent*)
{
  AnticoincidenceCollection = new GammaRayTelAnticoincidenceHitsCollection
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

G4bool GammaRayTelAnticoincidenceSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
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
  
  //  G4cout << ACDTileName << " " << edep/keV << G4endl;

  if (ACDTileName == "ACT" )
    // The hit is on an top ACD tile (ACDType 0)
    {
      // This is a new hit
      if (HitTopID[ACDTileNumber]==-1)
	{       
	  GammaRayTelAnticoincidenceHit* 
	    AnticoincidenceHit = new GammaRayTelAnticoincidenceHit;
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
	  GammaRayTelAnticoincidenceHit* 
	    AnticoincidenceHit = new GammaRayTelAnticoincidenceHit;
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
	  GammaRayTelAnticoincidenceHit* 
	    AnticoincidenceHit = new GammaRayTelAnticoincidenceHit;
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

void GammaRayTelAnticoincidenceSD::EndOfEvent(G4HCofThisEvent* HCE)
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

void GammaRayTelAnticoincidenceSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....














