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
// $Id: XrayFluoSD.cc
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
// 29 Nov 2002 Energy deposition bug fixed (Alfonso.mantero@ge.infn.it)
// 17 Jul 2003 Name changed to XrayFluoSD
// 01 Sep 2003 Constructor overload for different geometries handling
// -------------------------------------------------------------------

#include "XrayFluoSD.hh"
#include "XrayFluoSensorHit.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSD::XrayFluoSD(G4String name,
                                   XrayFluoDetectorConstruction* det)
  :G4VSensitiveDetector(name),Detector(0),planeDetector(0),mercuryDetector(0)
{

  Detector = det;
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[500];
  G4cout << "XrayFluoSD created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSD::XrayFluoSD(G4String name,
		       XrayFluoPlaneDetectorConstruction* det)
  :G4VSensitiveDetector(name),Detector(0),planeDetector(0),mercuryDetector(0)
{
  planeDetector = det;
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[500];
  G4cout << "XrayFluoSD created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSD::XrayFluoSD(G4String name,
		       XrayFluoMercuryDetectorConstruction* det)
  :G4VSensitiveDetector(name),Detector(0),planeDetector(0),mercuryDetector(0)
{
  mercuryDetector = det;
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[500];
  G4cout << "XrayFluoSD created" << G4endl;
}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




XrayFluoSD::~XrayFluoSD()
{

  delete [] HitHPGeID;

  // delete HPGeCollection;

  G4cout << "XrayFluoSD deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSD::Initialize(G4HCofThisEvent*)
 
//initializes HCE with the hits collection(s) created by this 
  //sensitive detector
{ 
  HPGeCollection = new XrayFluoSensorHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  
  G4int nPixel = 0;
  
  if (Detector) {nPixel = Detector->GetNbOfPixels();}
  else if (planeDetector) {nPixel = planeDetector->GetNbOfPixels();}
  else if (mercuryDetector) {nPixel = mercuryDetector->GetNbOfPixels();}
  
  for (G4int j=0;j<nPixel;j++)
    {HitHPGeID [j]= -1;};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XrayFluoSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{	



  G4double edep = aStep->GetTotalEnergyDeposit();

  /*
    G4String particleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();

    G4Track* track = aStep->GetTrack();
    G4int trackId = track->GetTrackID();
    //G4String processName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();

    G4double partEnergy = aStep->GetPreStepPoint()->GetKineticEnergy(); 
    G4double secondEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();
    G4cout << " la particella che deposita e': " << particleName << " ha una energia di keV "
	   << partEnergy  << " e deposita "<< edep << G4endl;
    G4cout << " la particella creata ha energia cinetica: " << secondEnergy << G4endl;  
  */

  if ((edep==0.)) return false;      

  //G4cout << " edep =  " << edep << G4endl;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    


  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  G4int PixelNumber = 0;

  if (Detector && Detector->GetNbOfPixels()>1) {PixelNumber= physVol->GetCopyNo();}
  else if (planeDetector && planeDetector->GetNbOfPixels()>1) {PixelNumber= physVol->GetCopyNo();}
  else if (mercuryDetector && mercuryDetector->GetNbOfPixels()>1) {PixelNumber= physVol->GetCopyNo();}


  
  if ( HitHPGeID[PixelNumber]==-1)
    { 
      XrayFluoSensorHit* HPGeHit = new XrayFluoSensorHit();
      HPGeHit->AddEnergy(edep);
      HitHPGeID[PixelNumber] = HPGeCollection->insert(HPGeHit) - 1;
      if (verboseLevel>0){
	G4cout << " New Hit on pixel: " << PixelNumber << G4endl;
      }
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

void XrayFluoSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,HPGeCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....















