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
//      ------------ GammaRayTelRunAction  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
// 18.11.2001 G.Santin
// - Modified the analysis management according to the new design
//
// ************************************************************

#include "GammaRayTelRunAction.hh"
#include "GammaRayTelAnalysis.hh"

#include <cstdlib>

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Threading.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelRunAction::GammaRayTelRunAction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelRunAction::~GammaRayTelRunAction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::BeginOfRunAction(const G4Run *run) {
    fRunID = run->GetRunID();

    // Master mode or sequential
    if (IsMaster()) {
        G4cout << "### Run " << run->GetRunID() << " starts (master)." << G4endl;
    } else {
        G4cout << "### Run " << run->GetRunID() << " starts (worker)." << G4endl;
    }

    // Prepare the visualization
    if (G4VVisManager::GetConcreteInstance() != nullptr) {
        auto *userInterface = G4UImanager::GetUIpointer();
        userInterface->ApplyCommand("/vis/scene/notifyHandlers");
    }

    // If analysis is used reset the histograms
    auto *analysis = GammaRayTelAnalysis::getInstance();
    // analysis->BeginOfRun(run->GetRunID());
    analysis->BeginOfRun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::EndOfRunAction(const G4Run *run) {
    G4cout << "End of Run " << run->GetRunID() << G4endl;

    // Close the file with the hits information
#ifdef G4STORE_DATA
	if (outFile != nullptr) {
		G4cout << "File " << fileName << G4endl;
		outFile->close();
		delete outFile;
		outFile = nullptr;
	}
#endif

    // If analysis is used, print out the histograms
    auto *analysis = GammaRayTelAnalysis::getInstance();
    analysis->EndOfRun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto GammaRayTelRunAction::GetOutputFile() -> std::ofstream* {
    if (outFile == nullptr) {
        OpenFile();
    }
    return outFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelRunAction::OpenFile() {
    // Open the file for the tracks of this run
#ifdef G4STORE_DATA
	// check that we are in a worker: returns -1 in a master and -2 in sequential
	// one file per thread is produced
	// Tracks_runR.N.dat, where R = run number, N = thread ID
    std::stringstream name;

    if (G4Threading::G4GetThreadId() >= 0) {
    	name << "Tracks_run" << fRunID << "." << G4Threading::G4GetThreadId() << ".dat";
    } else {
    	name << "Tracks_run" << fRunID << "." << G4Threading::G4GetThreadId() << ".dat";
    }

    if (outFile == nullptr) {
		outFile = new std::ofstream;
		outFile->open(name.str());
		fileName = G4String(name.str());
    }
    G4cout << "Open file: " << fileName << G4endl;
#endif
}
