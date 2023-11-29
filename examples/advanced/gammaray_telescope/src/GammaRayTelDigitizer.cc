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
//      ------------ GammaRayTelDigitizer  ------
//           by   F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include <vector>

#include "GammaRayTelAnticoincidenceHit.hh"
#include "GammaRayTelCalorimeterHit.hh"
#include "GammaRayTelDigi.hh"
#include "GammaRayTelDigitizer.hh"
#include "GammaRayTelDigitizerMessenger.hh"
#include "GammaRayTelTrackerHit.hh"

#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizer::GammaRayTelDigitizer(G4String name) : G4VDigitizerModule(name) {
    constexpr auto DIGIT_COLLECTION_NAME{"DigitsCollection"};
    collectionName.push_back(DIGIT_COLLECTION_NAME);

    // create a messenger for this class
    digitizerMessenger = new GammaRayTelDigitizerMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizer::~GammaRayTelDigitizer() {
    delete digitizerMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDigitizer::Digitize() {
    digitsCollection = new GammaRayTelDigitsCollection("GammaRayTelDigitizer", "DigitsCollection"); // to create the Digi Collection

    auto *digitManager = G4DigiManager::GetDMpointer();

    G4int trackerCollectionIdentifier{0}; // TrackerCollection
    G4int calorimeterCollectionIdentifier{0}; // CalorimeterCollection
    G4int anticoincidenceCollectionIdentifier{0}; // AnticoincidenceCollection

    totalEnergy = 0.; // at each event

    // TKR hits collection

    trackerCollectionIdentifier = digitManager->GetHitsCollectionID("TrackerCollection");
    GammaRayTelTrackerHitsCollection *trackerCollection{nullptr};
    trackerCollection = (GammaRayTelTrackerHitsCollection*) (digitManager->GetHitsCollection(trackerCollectionIdentifier));

    // CAL hits collection

    calorimeterCollectionIdentifier = digitManager->GetHitsCollectionID("CalorimeterCollection");
    GammaRayTelCalorimeterHitsCollection *calorimeterCollection{nullptr};
    calorimeterCollection = (GammaRayTelCalorimeterHitsCollection*) (digitManager->GetHitsCollection(calorimeterCollectionIdentifier));

    // ACD hits collection

    anticoincidenceCollectionIdentifier = digitManager->GetHitsCollectionID("AnticoincidenceCollection");
    GammaRayTelAnticoincidenceHitsCollection *anticoincidenceCollection{nullptr};
    anticoincidenceCollection = (GammaRayTelAnticoincidenceHitsCollection*) (digitManager->GetHitsCollection(anticoincidenceCollectionIdentifier));

    if (trackerCollection != nullptr) {
        G4int numberOfHits = trackerCollection->entries();

        for (auto i = 0; i < numberOfHits; i++) {
            G4double depositedEnergy = (*trackerCollection)[i]->GetDepositedEnergy();
            G4int stripNumber = (*trackerCollection)[i]->GetStripNumber();
            G4int planeNumber = (*trackerCollection)[i]->GetSiliconPlaneNumber();
            G4int IsX = (*trackerCollection)[i]->GetPlaneType();

            // digi generation only if energy deposit is greater than threshold

            if (depositedEnergy > energyThreshold) {
                auto *digit = new GammaRayTelDigi();
                digit->SetPlaneNumber(planeNumber);
                digit->SetPlaneType(IsX);
                digit->SetStripNumber(stripNumber);
                digit->SetDigitType(0);
                digit->SetEnergy(0.);
                digitsCollection->insert(digit);
            }
        }
    }

    if (calorimeterCollection != nullptr) {
        G4int numberOfHits = calorimeterCollection->entries();

        for (auto i = 0; i < numberOfHits; i++) {
            totalEnergy += (*calorimeterCollection)[i]->GetCALDepositedEnergy();
        }

        // digit generation only if energy deposit is greater than 0.

        if (totalEnergy > 0.) {
            auto *digit = new GammaRayTelDigi();
            digit->SetPlaneNumber(0);
            digit->SetPlaneType(0);
            digit->SetStripNumber(0);
            digit->SetDigitType(1);
            digit->SetEnergy(totalEnergy);
            digitsCollection->insert(digit);
        }
    }

    if (anticoincidenceCollection != nullptr) {
        G4int numberOfHits = anticoincidenceCollection->entries();

        for (G4int i = 0; i < numberOfHits; i++) {
            auto energy = (*anticoincidenceCollection)[i]->GetEdepACD();
            auto type = (*anticoincidenceCollection)[i]->GetACDTileNumber();

            // digit generation only if energy deposit is greater than 0.

            if (energy > acdThreshold) {
                auto *digit = new GammaRayTelDigi();
                digit->SetPlaneNumber(0);
                digit->SetPlaneType(0);
                digit->SetStripNumber(type);
                digit->SetDigitType(2);
                digit->SetEnergy(energy);
                digitsCollection->insert(digit);
            }
        }
    }

    if (trackerCollection != nullptr || anticoincidenceCollection != nullptr || calorimeterCollection != nullptr) {
        G4cout << "Number of digits in this event: " << digitsCollection->entries()
        // << " " << digitsCollection->GetName()
        // << " " << digitsCollection->GetDMname()
            << G4endl;
    }

    StoreDigiCollection(digitsCollection);

    G4int digitsCollectionIdentifier{-1};
    if (digitsCollectionIdentifier < 0) {
        // digiManager->List();
        digitsCollectionIdentifier = digitManager->GetDigiCollectionID("GammaRayTelDigitizer/DigitsCollection");
    }
}
