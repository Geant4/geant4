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


#include "FCALHadModuleSD.hh"

#include "FCALCalorHit.hh"

#include "FCALTestbeamSetup.hh"
#include "FCALHadModule.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALHadModuleSD::FCALHadModuleSD(G4String name) : G4VSensitiveDetector(name),
						  InitF2(0)
{
   HadModule = new FCALHadModule();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALHadModuleSD::~FCALHadModuleSD()
{
  delete HadModule;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::Initialize(G4HCofThisEvent*)
{
  if(InitF2 == 0) {
    HadModule->InitializeGeometry();
    InitF2++;
  }
  for (G4int j=0; j<2330; j++) { EvisF2Tile[j]=0.;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCALHadModuleSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep==0.) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());  
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();


    if(strcmp(physVol->GetName(),"F2LArGapPhysical")==0){
    G4int F2LArGapId = physVol->GetCopyNo();
    G4int F2TileId = HadModule->GetF2TileID(F2LArGapId);
    EvisF2Tile[F2TileId] = EvisF2Tile[F2TileId] + edep;
    };

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int NF2Tile = 0;
  G4int i=0;
  for (i=0; i<2330; i++){
    if(EvisF2Tile[i] > 0.) {
      NF2Tile++;
    };};

  G4cout << "Number of F2 tiles with Positive energy : " << NF2Tile <<  G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALHadModuleSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

