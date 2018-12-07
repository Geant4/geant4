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
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraPMTSD.cc
//    ****************************************************
//
//    Class used to define the Ultra photomultiplier as a sensitive detector.
//    Hits in this sensitive detector are defined in the UltraOpticalHit class
//
#include "UltraPMTSD.hh"

#include "G4Material.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraPMTSD::UltraPMTSD(G4String name):G4VSensitiveDetector(name)
{

  collectionName.insert("OpticalHitsCollection");

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraPMTSD::~UltraPMTSD()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraPMTSD::Initialize(G4HCofThisEvent* HCE)
{


  static int HCID1 = -1;


  // SensitiveDetectorName and collectionName are data members of G4VSensitiveDetector

  OpticalHitsCollection = 
    new UltraOpticalHitsCollection(SensitiveDetectorName,collectionName[0]);

  if(HCID1<0)
    { HCID1 = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID1,OpticalHitsCollection);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool UltraPMTSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{ 
 
  // Get Material 
  
  G4String thisVolume = aStep->GetTrack()->GetVolume()->GetName() ;
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();


  if (thisVolume != "PMT1" && thisVolume != "PMT2") 
    return false;
  if (particle != G4OpticalPhoton::Definition() ) 
    return false;

  if(particle == G4OpticalPhoton::Definition()) 
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);

  G4double      kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
  G4ThreeVector HitPosition   = aStep->GetPreStepPoint()->GetPosition() ;

  UltraOpticalHit* OpticalHit = new UltraOpticalHit ;
  OpticalHit->SetEnergy(kineticEnergy);
  OpticalHit->SetPosition(HitPosition);

 
  OpticalHitsCollection->insert(OpticalHit);


#ifdef ULTRA_VERBOSE
  G4cout << "*******************************" << G4endl;
  G4cout << "             PMT HIT           " << G4endl;
  G4cout << "  Volume:                      " << thisVolume << G4endl;
  G4cout << "  Photon energy (eV) :         " << kineticEnergy/CLHEP::eV << G4endl;
  G4cout << "  POSITION (mm) :              " 
	 << HitPosition.x()/CLHEP::mm << " " << HitPosition.y()/CLHEP::mm << " " << HitPosition.z()/CLHEP::mm << G4endl;
  G4cout << "*******************************" << G4endl;
#endif


  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraPMTSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection(HCID,OpticalHitsCollection);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraPMTSD::clear()
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraPMTSD::DrawAll()
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraPMTSD::PrintAll()
{;} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
