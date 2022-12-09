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

#include "GammaRayTelTrackerSD.hh"
#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerSD::GammaRayTelTrackerSD(G4String name) : G4VSensitiveDetector(name) {
    auto *runManager = G4RunManager::GetRunManager();

    detector = (GammaRayTelDetectorConstruction*) (runManager->GetUserDetectorConstruction());

    auto numberOfTKRTiles = detector->GetNbOfTKRTiles();
    numberOfTKRStrips = detector->GetNbOfTKRStrips();
    numberOfTKRLayers = detector->GetNbOfTKRLayers();
    numberOfTKRStrips = numberOfTKRStrips * numberOfTKRTiles;
    numberOfTKRChannels = numberOfTKRStrips * numberOfTKRTiles * numberOfTKRLayers;

    tkrHitXID = new G4int[numberOfTKRChannels];
    tkrHitYID = new G4int[numberOfTKRChannels];
    
    constexpr auto TRACKER_COLLECTION_NAME = "TrackerCollection";
    collectionName.insert(TRACKER_COLLECTION_NAME);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelTrackerSD::~GammaRayTelTrackerSD() {
    delete[] tkrHitXID;
    delete[] tkrHitYID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::Initialize(G4HCofThisEvent*) {
    trackerCollection = new GammaRayTelTrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

    for (auto i = 0; i < numberOfTKRChannels; i++) {
        tkrHitXID[i] = -1;
        tkrHitYID[i] = -1;
    };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelTrackerSD::ProcessHits(G4Step *step, G4TouchableHistory*) -> G4bool {
    G4double depositedEnergy = 0.;
    depositedEnergy = step->GetTotalEnergyDeposit();
    if (depositedEnergy == 0.) {
        return false;
    }

    auto totalNumberOfStrips = detector->GetNbOfTKRStrips();
    auto totalnumberOfTiles = detector->GetNbOfTKRTiles();

    // This TouchableHistory is used to obtain the physical volume of the hit
    auto *touchable = (G4TouchableHistory*) (step->GetPreStepPoint()->GetTouchable());
    auto *plane = touchable->GetVolume(2);

    G4int planeNumber = 0;
    planeNumber = plane->GetCopyNo();
    auto planeName = plane->GetName();

    // The hits sees now the real strip

    G4int stripNumber = 0;
    G4VPhysicalVolume *strip{nullptr};
    strip = touchable->GetVolume();

    G4String stripName = strip->GetName();
    stripNumber = strip->GetCopyNo();

    auto *tile = touchable->GetVolume(1);
    auto tileNumber = tile->GetCopyNo();
    auto tileName = tile->GetName();
    auto NTile = (tileNumber % totalnumberOfTiles);
    
    G4int channelNumber = 0;

    for (auto j = 0; j < totalnumberOfTiles; j++) {
        if (NTile == j) {
            stripNumber += totalNumberOfStrips * NTile;
        }
    }

    channelNumber = planeNumber * totalnumberOfTiles * totalNumberOfStrips + stripNumber;

/*
    G4cout << " Channel: " << channelNumber << G4endl;
    G4cout << " Plane: " << planeNumber << " " << planeName << G4endl;
    G4cout << " Strip: " << stripNumber << " " << stripName << G4endl;
*/

    // The hit is on an X silicon plane
    if (planeName == "TKRDetectorX") {
        if (tkrHitXID[channelNumber] == -1) { // This is a new hit
            auto *trackerHit = new GammaRayTelTrackerHit;
            trackerHit->SetPlaneType(1);
            trackerHit->AddDepositedEnergy(depositedEnergy);
            trackerHit->SetPosition(step->GetPreStepPoint()->GetPosition());
            trackerHit->SetSiliconPlaneNumber(planeNumber);
            trackerHit->SetStripNumber(stripNumber);
            tkrHitXID[channelNumber] = trackerCollection->insert(trackerHit) - 1;
        } else { // This is not new
            (*trackerCollection)[tkrHitXID[channelNumber]]->AddDepositedEnergy(depositedEnergy);
            // G4cout << "X" << planeNumber << " " << stripNumber << G4endl;
        }
    }

    // The hit is on an Y silicon plane
    if (planeName == "TKRDetectorY") {
        if (tkrHitYID[channelNumber] == -1) { // This is a new hit
            auto *trackerHit = new GammaRayTelTrackerHit;
            trackerHit->SetPlaneType(0);
            trackerHit->AddDepositedEnergy(depositedEnergy);
            trackerHit->SetPosition(step->GetPreStepPoint()->GetPosition());
            trackerHit->SetSiliconPlaneNumber(planeNumber);
            trackerHit->SetStripNumber(stripNumber);
            tkrHitYID[channelNumber] = trackerCollection->insert(trackerHit) - 1;
        } else { // This is not new
            (*trackerCollection)[tkrHitYID[channelNumber]]->AddDepositedEnergy(depositedEnergy);
            // G4cout << "Y" << planeNumber << " " << stripNumber << G4endl;
        }
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::EndOfEvent(G4HCofThisEvent *collection) {
    static G4int collectionIdentifier = -1;
    if (collectionIdentifier < 0) {
        collectionIdentifier = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    collection->AddHitsCollection(collectionIdentifier, trackerCollection);

    for (auto i = 0; i < numberOfTKRChannels; i++) {
        tkrHitXID[i] = -1;
        tkrHitYID[i] = -1;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::clear() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::DrawAll() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelTrackerSD::PrintAll() {
}
