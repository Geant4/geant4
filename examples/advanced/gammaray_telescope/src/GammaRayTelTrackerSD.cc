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
//      ------------ GammaRayTelTrackerSD  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
#include "G4RunManager.hh"
#include "GammaRayTelTrackerSD.hh"

#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerSD::GammaRayTelTrackerSD(G4String name):G4VSensitiveDetector(name)
{
 G4RunManager* runManager = G4RunManager::GetRunManager();
  Detector =
    (GammaRayTelDetectorConstruction*)(runManager->GetUserDetectorConstruction());
  
  G4int NbOfTKRTiles  =  Detector->GetNbOfTKRTiles();
  NbOfTKRStrips  = Detector->GetNbOfTKRStrips();
  NbOfTKRLayers  = Detector->GetNbOfTKRLayers();  
  NbOfTKRStrips = NbOfTKRStrips*NbOfTKRTiles;  
  
  NbOfTKRChannels = NbOfTKRStrips* NbOfTKRTiles * NbOfTKRLayers;
  
  ThitXID = new G4int[NbOfTKRChannels];
  ThitYID = new G4int[NbOfTKRChannels];
  collectionName.insert("TrackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerSD::~GammaRayTelTrackerSD()
{
  delete [] ThitXID;
  delete [] ThitYID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::Initialize(G4HCofThisEvent*)
{
  TrackerCollection = new GammaRayTelTrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);

 for (G4int i=0;i<NbOfTKRChannels;i++)
   {
     ThitXID[i] = -1;
     ThitYID[i] = -1;
   };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool GammaRayTelTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* )
{ 
   
  G4double edep = 0.;
  edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;      
  
  G4int StripTotal = Detector->GetNbOfTKRStrips();
  G4int TileTotal  = Detector->GetNbOfTKRTiles();  

  // This TouchableHistory is used to obtain the physical volume
  // of the hit
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  
  //G4VPhysicalVolume* phys_tile = theTouchable->GetVolume();  
  
  G4VPhysicalVolume* plane = theTouchable->GetVolume(2);  
  
  G4int PlaneNumber = 0;
  PlaneNumber=plane->GetCopyNo();
  G4String PlaneName = plane->GetName();

  // The hits sees now the real strip

  G4int StripNumber = 0;
  G4VPhysicalVolume* strip = 0;
  strip = theTouchable->GetVolume();

  G4String StripName = strip->GetName();
  StripNumber= strip->GetCopyNo();  
  
  G4VPhysicalVolume* tile = theTouchable->GetVolume(1); 
  G4int TileNumber = tile->GetCopyNo();  
  G4String TileName = tile->GetName();   
  
  G4int NTile = (TileNumber%TileTotal);  
  G4int j=0;
  
  G4int NChannel = 0;
  
  for (j=0;j<TileTotal;j++)
    {
      if(NTile==j) StripNumber += StripTotal*NTile;
    }  
  
  NChannel = PlaneNumber*TileTotal*StripTotal + StripNumber;

  /*  G4cout << NChannel << " Channel Number" << G4endl;
      G4cout << " Plane Number = " << PlaneNumber << " " << PlaneName 
      << G4endl;
      G4cout << StripName << " " << StripNumber << G4endl;      */ 
  
  if (PlaneName == "TKRDetectorX" )
    // The hit is on an X silicon plane
    {
      // This is a new hit
      if (ThitXID[NChannel]==-1)
	{       
	  GammaRayTelTrackerHit* TrackerHit = new GammaRayTelTrackerHit;
	  TrackerHit->SetPlaneType(1);
	  TrackerHit->AddSil(edep);
	  TrackerHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  TrackerHit->SetNSilPlane(PlaneNumber);
	  TrackerHit->SetNStrip(StripNumber);
	  ThitXID[NChannel] = 
	    TrackerCollection->insert(TrackerHit) -1;
	}
      else // This is not new
	{
	  (*TrackerCollection)[ThitXID[NChannel]]->AddSil(edep);
          // G4cout << "X" << PlaneNumber << " " << StripNumber << G4endl;
	}
    }
  
  if (PlaneName == "TKRDetectorY")
    // The hit is on an Y silicon plane    
    {   
      // This is a new hit
      if (ThitYID[NChannel]==-1)
	{       
	  GammaRayTelTrackerHit* TrackerHit = new GammaRayTelTrackerHit;
	  TrackerHit->SetPlaneType(0);
	  TrackerHit->AddSil(edep);
	  TrackerHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	  TrackerHit->SetNSilPlane(PlaneNumber);
	  TrackerHit->SetNStrip(StripNumber);
	  ThitYID[NChannel] = 
	    TrackerCollection->insert(TrackerHit)-1;
	}
      else // This is not new
	{
	  (*TrackerCollection)[ThitYID[NChannel]]->AddSil(edep);
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


  for (G4int i=0;i<NbOfTKRChannels;i++) 
    {
      ThitXID[i] = -1;
      ThitYID[i] = -1;
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













