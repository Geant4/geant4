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
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
// 29 Nov 2002 Energy deposition bug fixed (Alfonso.mantero@ge.infn.it)
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
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoHPGeSD::XrayFluoHPGeSD(G4String name,
                                   XrayFluoDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det)
{
  collectionName.insert("HPGeCollection");
  HitHPGeID = new G4int[500];
  //G4cout << "XrayFluoHPGeSD created" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoHPGeSD::~XrayFluoHPGeSD()
{

  delete [] HitHPGeID;

  // delete HPGeCollection;

  G4cout << "XrayFluoHPGeSD deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoHPGeSD::Initialize(G4HCofThisEvent*HCE)
 
//initializes HCE with the hits collection(s) created by this 
  //sensitive detector
{ 
  if (HCE){
    
    HPGeCollection = new XrayFluoSensorHitsCollection
      (SensitiveDetectorName,collectionName[0]); 
    for (G4int j=0;j<Detector->GetNbOfPixels();j++)
      {HitHPGeID [j]= -1;};
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool XrayFluoHPGeSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{	
 
  if (ROhist){
    
    G4double edep = aStep->GetTotalEnergyDeposit();
    //    G4String particleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
    //    G4String processName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
    
    //    G4double partEnergy = aStep->GetPreStepPoint()->GetKineticEnergy(); 
    //    G4double secondEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();
    //    G4cout << " la particella che deposita e': " << particleName << " ha una energia di keV "
    //  	 << partEnergy  << " e deposita "<< edep << G4endl;
    //    G4cout << " la particella creata deposita: " << secondEnergy << G4endl;  
    //    if (processName){
    
    
    
    //      G4cout << " la particella viene da: " << processName << G4endl;
    //    }
    
    //    else {G4cout << " il deposito proviene da: Rayleugh" << G4endl;}
    
    if (edep==0.) return false;      
    
    //  G4cout << " edep =  " << edep << G4endl;
    
    G4TouchableHistory* theTouchable
      = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
    G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
    //  theTouchable->MoveUpHistory();     
    G4int PixelNumber = 0;
    if (Detector->GetNbOfPixels()>1) PixelNumber= physVol->GetCopyNo() ;
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
  else {return false;}

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















