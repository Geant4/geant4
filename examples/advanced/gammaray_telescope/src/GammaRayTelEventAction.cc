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
//      ------------ GammaRayTelEventAction  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// - inclusion of Digits by F.Longo & R.Giannitrapani (24 oct 2001)
// 
// - Modification of analysis management by G.Santin (18 Nov 2001)
// 
// ************************************************************

#include "GammaRayTelEventAction.hh"
#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelAnticoincidenceHit.hh"
#include "GammaRayTelCalorimeterHit.hh"
#include "GammaRayTelAnalysis.hh"
#include "GammaRayTelDigi.hh"
#include "GammaRayTelDigitizer.hh"

#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VHitsCollection.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This file is a global variable in which we store energy deposition per hit
// and other relevant information

GammaRayTelEventAction::GammaRayTelEventAction(GammaRayTelRunAction *runAction) : theRunAction(runAction) {
    auto *digitizer = new GammaRayTelDigitizer("GammaRayTelDigitizer");
    G4DigiManager::GetDMpointer()->AddNewModule(digitizer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::~GammaRayTelEventAction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::BeginOfEventAction(const G4Event *event) {
    G4int eventIdentifer = event->GetEventID();
    G4cout << "Event: " << eventIdentifer << G4endl;
    auto *sensitiveDetectorManager = G4SDManager::GetSDMpointer();

    if (trackerCollectionID == -1) {
        trackerCollectionID = sensitiveDetectorManager->GetCollectionID("TrackerCollection");
    }
    if (anticoincidenceCollectionID == -1) {
        anticoincidenceCollectionID = sensitiveDetectorManager->GetCollectionID("AnticoincidenceCollection");
    }
    if (calorimeterCollectionID == -1) {
        calorimeterCollectionID = sensitiveDetectorManager->GetCollectionID("CalorimeterCollection");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::EndOfEventAction(const G4Event *event) {
    G4int eventIdentifier = event->GetEventID();

    if (theRunAction == nullptr) {
        G4Exception("GammaRayTelEventAction::BeginOfEventAction()", "GTR0001", FatalException, "Null pointer to Run Action: this should not be");
    }

#ifdef G4STORE_DATA
    auto *outputFile = theRunAction->GetOutputFile();
#endif

    auto *collection = event->GetHCofThisEvent();
    GammaRayTelTrackerHitsCollection *trackerCollection{nullptr}; // TKR
    // GammaRayTelCalorimeterHitsCollection* calorimeterCollection{nullptr}; // CAL
    // GammaRayTelAnticoincidenceHitsCollection* anticoincidenceCollection{nullptr}; // ACD

    auto *fDM = G4DigiManager::GetDMpointer();
    GammaRayTelAnalysis *analysis = GammaRayTelAnalysis::getInstance();

    if (collection != nullptr) {
        trackerCollection = (GammaRayTelTrackerHitsCollection*) (collection->GetHC(trackerCollectionID));
        // calorimeterCollection = (GammaRayTelCalorimeterHitsCollection*) (collection->GetHC(calorimeterCollectionID));
        // anticoincidenceCollection = (GammaRayTelAnticoincidenceHitsCollection*) (collection->GetHC(anticoincidenceCollectionID));

        if (trackerCollection != nullptr) {
            G4int numberOfHits = trackerCollection->entries();
            G4cout << "Number of tracker hits in this event: " << numberOfHits << G4endl;

            G4double depositedEnergy{0};
            G4int stripNumber;
            G4int planeNumber;
            G4int isXPlane;

            // This is a cycle on all the tracker hits of this event

            for (auto i = 0; i < numberOfHits; i++) {
                // Here we put the hit data in an ASCII file for later analysis
                depositedEnergy = (*trackerCollection)[i]->GetDepositedEnergy();
                stripNumber = (*trackerCollection)[i]->GetStripNumber();
                planeNumber = (*trackerCollection)[i]->GetSiliconPlaneNumber();
                isXPlane = (*trackerCollection)[i]->GetPlaneType();

#ifdef G4STORE_DATA
                (*outputFile) << std::setw(7) << eventIdentifier
                    << " " << depositedEnergy/keV
                    << " " << stripNumber
                    << " " << planeNumber
                    << " " << isXPlane
                    << " " << (*trackerCollection)[i]->GetPosition().x() / mm
                    << " " << (*trackerCollection)[i]->GetPosition().y() / mm
                    << " " << (*trackerCollection)[i]->GetPosition().z() / mm
                    << " " << G4endl;
#else 	  
                G4cout << std::setw(7) << eventIdentifier
                    << " " << depositedEnergy / keV
                    << " " << stripNumber
                    << " " << planeNumber
                    << " " << isXPlane
                    << " " << (*trackerCollection)[i]->GetPosition().x() / mm
                    << " " << (*trackerCollection)[i]->GetPosition().y() / mm
                    << " " << (*trackerCollection)[i]->GetPosition().z() / mm
                    << " " << G4endl;
#endif

                // Here we fill the histograms of the Analysis manager
                if (isXPlane != 0) {
                    if (analysis->GetHisto2DMode() == "position") {
                        analysis->InsertPositionXZ((*trackerCollection)[i]->GetPosition().x() / mm, (*trackerCollection)[i]->GetPosition().z() / mm);
                    } else {
                        analysis->InsertPositionXZ(stripNumber, planeNumber);
                    }

                    if (planeNumber == 0) {
                        analysis->InsertEnergy(depositedEnergy / keV);
                    }
                    analysis->InsertHits(planeNumber);
                } else {
                    if (analysis->GetHisto2DMode() == "position") {
                        analysis->InsertPositionYZ((*trackerCollection)[i]->GetPosition().y() / mm, (*trackerCollection)[i]->GetPosition().z() / mm);
                    } else {
                        analysis->InsertPositionYZ(stripNumber, planeNumber);
                    }
                    if (planeNumber == 0) {
                        analysis->InsertEnergy(depositedEnergy / keV);
                    }
                    analysis->InsertHits(planeNumber);
                }
                analysis->setNtuple(depositedEnergy / keV, planeNumber,
                    (*trackerCollection)[i]->GetPosition().x() / mm,
                    (*trackerCollection)[i]->GetPosition().y() / mm,
                    (*trackerCollection)[i]->GetPosition().z() / mm
                );
            }
            analysis->EndOfEvent(numberOfHits);
        }

        auto *myDM = (GammaRayTelDigitizer*) fDM->FindDigitizerModule("GammaRayTelDigitizer");
        myDM->Digitize();

#ifdef G4STORE_DATA
        // The whole block is needed only when outputFile is active;
        // protect block to avoid compilations warnings from gcc4.6, Gunter Folger

		auto digitsCollectionIdentifier = fDM->GetDigiCollectionID("DigitsCollection");
        // G4cout << "Digits collection: " << digitsCollectionIdentifier << G4endl;
      
		auto *digitsCollection = (GammaRayTelDigitsCollection*) fDM->GetDigiCollection(digitsCollectionIdentifier);
      
        if (digitsCollection != nullptr) {
            auto numberOfDigits = digitsCollection->entries();
            // G4cout << "Total number of digits: " << numberOfDigits << G4endl;

            G4int stripNumber;
            G4int planeNumber;
            G4int isXPlane;

            for (auto i = 0; i < numberOfDigits; i++) {
                // Here we put the digi data in an ASCII file for later analysis
                stripNumber = (*digitsCollection)[i]->GetStripNumber();
                planeNumber = (*digitsCollection)[i]->GetPlaneNumber();
                isXPlane = (*digitsCollection)[i]->GetPlaneType();

                (*outputFile) << std::setw(7)
                    << eventIdentifier
                    << " " << stripNumber
                    << " " << planeNumber
                    << " " << isXPlane
                    << " " << G4endl;
            }
        }
#endif      
    }
}
