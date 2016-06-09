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

#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

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
G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();


if (thisVolume != "PMT1" && thisVolume != "PMT2") return false;
if (particleName != "opticalphoton" ) return false;

if(particleName == "opticalphoton") aStep->GetTrack()->SetTrackStatus(fStopAndKill);

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
 G4cout << "  Photon energy (eV) :         " << kineticEnergy/eV << G4endl;
 G4cout << "  POSITION (mm) :              " 
        << HitPosition.x()/mm << " " << HitPosition.y()/mm << " " << HitPosition.z()/mm << G4endl;
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
