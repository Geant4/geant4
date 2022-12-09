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
//      ------------ GammaRayAnalysisManager  ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
//
// 03.04.2013 F.Longo/L.Pandola
// - migrated to G4tools
//
// 29.05.2003 F.Longo 
// - Anaphe 5.0.5 compliant
//
// 18.06.2002 R.Giannitrapani, F.Longo & G.Santin
// - new release for Anaphe 4.0.3
//
// 07.12.2001 A.Pfeiffer
// - integrated Guy's addition of the ntuple
//
// 06.12.2001 A.Pfeiffer
// - updating to new design (singleton)
//
// 22.11.2001 G.Barrand
// - Adaptation to AIDA
//
// ************************************************************

#include <fstream>
#include <iomanip>

#include "GammaRayTelAnalysis.hh"
#include "GammaRayTelAnalysisMessenger.hh"
#include "GammaRayTelDetectorConstruction.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis *GammaRayTelAnalysis::instance{nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis::GammaRayTelAnalysis() : detector(nullptr), histo2DMode("strip") {
    detector = static_cast<const GammaRayTelDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    // Define the messenger and the analysis system
    analysisMessenger = new GammaRayTelAnalysisMessenger(this);
    histogramFileName = "gammaraytel";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysis::~GammaRayTelAnalysis() {
    Finish();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysis::Init() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysis::Finish() {
    delete analysisMessenger;
    analysisMessenger = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelAnalysis::getInstance() -> GammaRayTelAnalysis* {
	if (instance == nullptr) {
		instance = new GammaRayTelAnalysis();
	}
	return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This function fill the 2d histogram of the XZ positions
void GammaRayTelAnalysis::InsertPositionXZ(G4double x, G4double z) {
	auto *manager = G4AnalysisManager::Instance();
	manager->FillH2(1, x, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This function fill the 2d histogram of the YZ positions
void GammaRayTelAnalysis::InsertPositionYZ(G4double y, G4double z) {
	auto *manager = G4AnalysisManager::Instance();
	manager->FillH2(2, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This function fill the 1d histogram of the energy released in the last Si plane
void GammaRayTelAnalysis::InsertEnergy(G4double energy) {
	auto *manager = G4AnalysisManager::Instance();
	manager->FillH1(1, energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This function fill the 1d histogram of the hits distribution along the TKR planes
void GammaRayTelAnalysis::InsertHits(G4int planeNumber) {
	auto *manager = G4AnalysisManager::Instance();
	manager->FillH1(2, planeNumber);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysis::setNtuple(G4double energy, G4int planeNumber, G4double x, G4double y, G4double z) {
	auto *manager = G4AnalysisManager::Instance();
	manager->FillNtupleDColumn(0, energy);
	manager->FillNtupleDColumn(1, planeNumber);
	manager->FillNtupleDColumn(2, x);
	manager->FillNtupleDColumn(3, y);
	manager->FillNtupleDColumn(4, z);
	manager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/* 
 This member reset the histograms and it is called at the begin of each run;
 here we put the inizialization so that the histograms have
 always the right dimensions depending from the detector geometry
 */
void GammaRayTelAnalysis::BeginOfRun() {
	auto *manager = G4AnalysisManager::Instance();
	manager->SetDefaultFileType("root");

	// Open an output file

	G4cout << "Opening output file " << histogramFileName << " ... ";
	manager->OpenFile(histogramFileName);
	manager->SetFirstHistoId(1);
	G4cout << " done" << G4endl;

	auto Nplane = detector->GetNbOfTKRLayers();
	auto numberOfStrips = detector->GetNbOfTKRStrips();
	auto numberOfTiles = detector->GetNbOfTKRTiles();
	auto sizeXY = detector->GetTKRSizeXY();
	auto sizeZ = detector->GetTKRSizeZ();
	auto N = numberOfStrips * numberOfTiles;

	// Book 1D histograms
	//-------------------

	constexpr auto NUMBER_OF_BINS{100};
	constexpr auto LOWER_BOUND{50};
	constexpr auto UPPER_BOUND{200};

	// 1D histogram that store the energy deposition of the particle in the last (number 0) TKR X-plane
	manager->CreateH1("1", "Deposited energy in the last X plane (keV)", NUMBER_OF_BINS, LOWER_BOUND, UPPER_BOUND);

	// 1D histogram that store the hits distribution along the TKR X-planes
	manager->CreateH1("2", "Hits distribution in TKR X planes", Nplane, 0, Nplane - 1);

	// Book 2D histograms 
	//-------------------

	// 2D histogram that store the position (mm) of the hits (XZ projection)

    if (histo2DMode == "strip") {
        manager->CreateH2("1", "Tracker hits XZ (strip, plane)", N, 0, N - 1, 2 * Nplane, 0, Nplane - 1);
    } else {
        manager->CreateH2("1", "Tracker hits XZ (x, z) in mm", G4int(sizeXY / 5), -sizeXY / 2, sizeXY / 2, G4int(sizeZ / 5), -sizeZ / 2, sizeZ / 2);
    }

    // 2D histogram that store the position (mm) of the hits (YZ projection)
    if (histo2DMode == "strip") {
        manager->CreateH2("2", "Tracker hits YZ (strip, plane)", N, 0, N - 1, 2 * Nplane, 0, Nplane - 1);
    } else {
        manager->CreateH2("2", "Tracker hits YZ (x, z) in mm", G4int(sizeXY / 5), -sizeXY / 2, sizeXY / 2, G4int(sizeZ / 5), -sizeZ / 2, sizeZ / 2);
    }

	// Book n-tuple (energy, plane, x, y, z)
	//------------------------------------------  
	manager->CreateNtuple("1", "Track n-tuple");
	manager->CreateNtupleDColumn("energy");
	manager->CreateNtupleDColumn("plane");
	manager->CreateNtupleDColumn("x");
	manager->CreateNtupleDColumn("y");
	manager->CreateNtupleDColumn("z");
	manager->FinishNtuple();
}

/* 
 This member is called at the end of each run 
 */
void GammaRayTelAnalysis::EndOfRun() {
	// Save histograms
	auto *manager = G4AnalysisManager::Instance();
	manager->Write();
	manager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This member is called at the end of every event
void GammaRayTelAnalysis::EndOfEvent(G4int /* flag */) {
}
