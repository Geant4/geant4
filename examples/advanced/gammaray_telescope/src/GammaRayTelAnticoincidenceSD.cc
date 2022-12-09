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

#include "GammaRayTelAnticoincidenceSD.hh"
#include "GammaRayTelAnticoincidenceHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceSD::GammaRayTelAnticoincidenceSD(G4String name) : G4VSensitiveDetector(name) {
    auto *runManager = G4RunManager::GetRunManager();
    detector = (GammaRayTelDetectorConstruction*) (runManager->GetUserDetectorConstruction());

    numberOfACDLateralTiles = detector->GetNbOfACDLateralTiles();
    numberOfACDTopTiles = detector->GetNbOfACDTopTiles();

    G4cout << numberOfACDLateralTiles << " LATERAL " << G4endl;
    G4cout << numberOfACDTopTiles << " TOP " << G4endl;

    hitLateralID = new G4int[numberOfACDLateralTiles];
    hitTopID = new G4int[numberOfACDTopTiles];
    collectionName.insert("AnticoincidenceCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnticoincidenceSD::~GammaRayTelAnticoincidenceSD() {
    delete[] hitLateralID;
    delete[] hitTopID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::Initialize(G4HCofThisEvent*) {
    anticoincidenceCollection = new GammaRayTelAnticoincidenceHitsCollection(SensitiveDetectorName, collectionName[0]);
    for (auto i = 0; i < numberOfACDLateralTiles; i++) {
        hitLateralID[i] = -1;
    }

    for (auto j = 0; j < numberOfACDTopTiles; j++) {
        hitTopID[j] = -1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelAnticoincidenceSD::ProcessHits(G4Step *step, G4TouchableHistory*) -> G4bool {
    auto depositedEnergy = step->GetTotalEnergyDeposit();
    if (depositedEnergy / keV == 0.) {
        return false;
    }

    // This TouchableHistory is used to obtain the physical volume of the hit
    auto *theTouchable = (G4TouchableHistory*) (step->GetPreStepPoint()->GetTouchable());

    auto *acdTile = theTouchable->GetVolume();
    auto acdTileNumber = acdTile->GetCopyNo();
    auto acdTileName = acdTile->GetName();

    // G4cout << acdTileName << " " << depositedEnergy / keV << G4endl;

    if (acdTileName == "ACT") { // The hit is on a top ACD tile (ACDType: 0)
        if (hitTopID[acdTileNumber] == -1) { // This is a new hit
            auto *anticoincidenceHit = new GammaRayTelAnticoincidenceHit;
            anticoincidenceHit->SetACDType(0);
            anticoincidenceHit->AddEnergy(depositedEnergy);
            anticoincidenceHit->SetPosition(step->GetPreStepPoint()->GetPosition());
            anticoincidenceHit->SetACDTileNumber(acdTileNumber);
            hitTopID[acdTileNumber] = anticoincidenceCollection->insert(anticoincidenceHit) - 1;
        } else { // This is not new
            (*anticoincidenceCollection)[hitTopID[acdTileNumber]]->AddEnergy(depositedEnergy);
        }
    }

    if (acdTileName == "ACL1") { // The hit is on a lateral (left-right) ACD tile (ACDType: 1)
        if (hitLateralID[acdTileNumber] == -1) { // This is a new hit
            auto *anticoincidenceHit = new GammaRayTelAnticoincidenceHit;
            anticoincidenceHit->SetACDType(1);
            anticoincidenceHit->AddEnergy(depositedEnergy);
            anticoincidenceHit->SetPosition(step->GetPreStepPoint()->GetPosition());
            anticoincidenceHit->SetACDTileNumber(acdTileNumber);
            hitLateralID[acdTileNumber] = anticoincidenceCollection->insert(anticoincidenceHit) - 1;
        } else { // This is not new
            (*anticoincidenceCollection)[hitLateralID[acdTileNumber]]->AddEnergy(depositedEnergy);
        }
    }

    if (acdTileName == "ACL2") { // The hit is on a lateral (rear-front) ACD tile (ACDType: 2)
        if (hitLateralID[acdTileNumber] == -1) { // This is a new hit
            auto *anticoincidenceHit = new GammaRayTelAnticoincidenceHit;
            anticoincidenceHit->SetACDType(2);
            anticoincidenceHit->AddEnergy(depositedEnergy);
            anticoincidenceHit->SetPosition(step->GetPreStepPoint()->GetPosition());
            anticoincidenceHit->SetACDTileNumber(acdTileNumber);
            hitLateralID[acdTileNumber] = anticoincidenceCollection->insert(anticoincidenceHit) - 1;
        } else { // This is not new
            (*anticoincidenceCollection)[hitLateralID[acdTileNumber]]->AddEnergy(depositedEnergy);
        }
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::EndOfEvent(G4HCofThisEvent *collection) {
    static G4int collectionIdentifier = -1;
    if (collectionIdentifier < 0) {
        collectionIdentifier = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    collection->AddHitsCollection(collectionIdentifier, anticoincidenceCollection);

    for (auto i = 0; i < numberOfACDLateralTiles; i++) {
        hitLateralID[i] = -1;
    }

    for (auto j = 0; j < numberOfACDTopTiles; j++) {
        hitTopID[j] = -1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::clear() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::DrawAll() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnticoincidenceSD::PrintAll() {
}
