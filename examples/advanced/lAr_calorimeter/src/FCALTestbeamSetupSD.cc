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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <iostream>

#include "FCALTestbeamSetupSD.hh"

#include "FCALCalorHit.hh"

#include "FCALTestbeamSetup.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTestbeamSetupSD::FCALTestbeamSetupSD(G4String name) : G4VSensitiveDetector(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALTestbeamSetupSD::~FCALTestbeamSetupSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::Initialize(G4HCofThisEvent*)
{
  EBeamS1 = EBeamS2 = EBeamS3 = 0.;
  EHoleScint = EBeamHole = 0.;
  EBeamDead = 0;
  G4int j;
  for (j =0 ; j<NLENGTH ; j++) { 
    ETailVis[j] = 0.;
    ETailDep[j] = 0.;
  }
  TailCatcherID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCALTestbeamSetupSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep==0.) return true;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4String name = physVol->GetName();
  TailCatcherID = physVol->GetCopyNo();
  
  if(name == "ScintS1Physical") { 
    EBeamS1 = EBeamS1 + edep;}
  
  else if(name == "ScintS2Physical") { 
    EBeamS2 = EBeamS2 + edep;}
  
  else if(name == "ScintS3Physical") {  
    EBeamS3 = EBeamS3 + edep;}

  else if(name == "HoleScintPhysical"){ EHoleScint += edep;}
  else if(name == "HoleCntrScintPhysical"){
    EBeamHole = EBeamHole + edep;}

  else if(name == "MWPCPhysical") { EBeamDead += edep;}
  else if(name == "HoleCntrPbPhysical") { EBeamDead += edep;}
  else if(name == "HoleCntrAlPhysical") { EBeamDead += edep;}
  else if(name == "LeadWallPhysical") { EBeamDead += edep;}
  else if(name == "IronWallPhysical") { EBeamDead += edep;}
  else if(TailCatcherID >= 0 && TailCatcherID < NLENGTH) {
    if(name == "BigScintPhysical") {
      ETailVis[TailCatcherID] += edep;
    }
    else if(name == "SmallScintPhysical") {
      if(TailCatcherID+3 < NLENGTH)  ETailVis[TailCatcherID + 3] += edep;
    }
    else if(name == "BigIronPhysical") {
      ETailDep[TailCatcherID] += edep;
    }
    else if(name == "SmallIronPhysical") {
      if(TailCatcherID+2 < NLENGTH)  ETailDep[TailCatcherID+2] += edep;
    }
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::EndOfEvent(G4HCofThisEvent*)
{
  G4cout << " Visisble Energy in S1 , S2 , S3 in (MeV)" << G4endl;
  G4cout << EBeamS1/MeV << " " << EBeamS2/MeV << " " << EBeamS3/MeV << " " << G4endl;

  G4cout << " Visible Energy in Hole Counter  (MeV) " << G4endl;
  G4cout << EHoleScint/MeV << " " << EBeamHole/MeV << G4endl;

  G4cout << " Visible Energy in Upstream Dead Materials " << G4endl;
  G4cout << EBeamDead/MeV << G4endl;

  G4cout << " Visible Energy in Tail Catcher Scintillator" << G4endl;
  G4int j;
  for (j=1; j<8 ; j++) {G4cout <<  ETailVis[j]/MeV << " " ;};  G4cout << G4endl;
 
  G4cout << " Visible Energy in Tail Catcher Absorber" << G4endl;
  for (j=1; j<7 ; j++) {G4cout <<  ETailDep[j]/MeV << " " ;};  G4cout << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALTestbeamSetupSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

