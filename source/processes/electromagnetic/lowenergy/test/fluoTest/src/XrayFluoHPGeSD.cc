//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoHPGeSD.cc
// GEANT4 tag $Name: 
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoHPGeSD.hh"
#include "XrayFluoSensorHit.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoHPGeSD::XrayFluoHPGeSD(G4String name,
                                   XrayFluoDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[500];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoHPGeSD::~XrayFluoHPGeSD()
{
  delete [] HitHPGeID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoHPGeSD::Initialize(G4HCofThisEvent*HCE)
 
//initializes HCE with the hits collection(s) created by this 
  //sensitive detector
{ 
  HPGeCollection = new XrayFluoSensorHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
 for (G4int j=0;j<Detector->GetNbOfPixels();j++)
 {HitHPGeID [j]= -1;};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XrayFluoHPGeSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
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
      XrayFluoSensorHit* HPGeHit = new XrayFluoSensorHit();
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

void XrayFluoHPGeSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,HPGeCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoHPGeSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoHPGeSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoHPGeSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




