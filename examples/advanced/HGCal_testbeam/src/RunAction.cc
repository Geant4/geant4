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
#include "RunAction.hh"
#include "EventAction.hh"

#include "G4UserRunAction.hh"
#include "G4GenericMessenger.hh"
#include "G4String.hh"
#include "G4AnalysisManager.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction *eventAction)
    : G4UserRunAction(), fEventAction(eventAction),
      fOutputFileDir("sim_HGCalTB_G4Standalone") {

  fMessenger =
      new G4GenericMessenger(this, "/HGCalTestbeam/output/", "Output control");

  // randomizePrimary command
  auto &fileNameCommand = fMessenger->DeclareProperty("file", fOutputFileDir);
  G4String guidance = "Define output file location.";
  fileNameCommand.SetGuidance(guidance);
  fileNameCommand.SetParameterName("filename", true);
  fileNameCommand.SetDefaultValue("sim_HGCalTB_G4Standalone");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *) {
  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespaces
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(1);
  std::cout << "Output file is: " << fOutputFileDir << std::endl;
  analysisManager->SetFileName(fOutputFileDir);

  if (fEventAction) {
    analysisManager->CreateNtuple("hits", "hits");
    analysisManager->CreateNtupleIColumn("event"); // column Id = 0
    analysisManager->CreateNtupleIColumn(
        "pdgID", fEventAction->fPrimariesPDG); // column Id = 1
    analysisManager->CreateNtupleDColumn(
        "beamEnergy", fEventAction->fPrimariesEnergy); // column Id = 2
    analysisManager->CreateNtupleDColumn(
        "beamX_cm", fEventAction->fPrimariesX); // column Id = 3
    analysisManager->CreateNtupleDColumn(
        "beamY_cm", fEventAction->fPrimariesY); // column Id = 4
    analysisManager->CreateNtupleDColumn(
        "beamZ_cm", fEventAction->fPrimariesZ); // column Id = 5
    analysisManager->CreateNtupleIColumn("siliconHits_ID",
                                         fEventAction->fSiHitsID);
    analysisManager->CreateNtupleDColumn("siliconHits_x_cm",
                                         fEventAction->fSiHitsX);
    analysisManager->CreateNtupleDColumn("siliconHits_y_cm",
                                         fEventAction->fSiHitsY);
    analysisManager->CreateNtupleDColumn("siliconHits_z_cm",
                                         fEventAction->fSiHitsZ);
    analysisManager->CreateNtupleDColumn("siliconHits_Edep_keV",
                                         fEventAction->fSiHitsEdep);
    analysisManager->CreateNtupleDColumn("siliconHits_EdepNonIonizing_keV",
                                         fEventAction->fSiHitsEdepNonIonising);
    analysisManager->CreateNtupleDColumn("siliconHits_TOA_ns",
                                         fEventAction->fSiHitsTOA);
    analysisManager->CreateNtupleDColumn("siliconHits_TOA_last_ns",
                                         fEventAction->fSiHitsTOA);
    analysisManager->CreateNtupleIColumn("siliconHits_type",
                                         fEventAction->fSiHitsType);

    analysisManager->CreateNtupleIColumn("SiPMHits_ID",
                                         fEventAction->fSiPMhitsID);
    analysisManager->CreateNtupleDColumn("SiPMHits_x_cm",
                                         fEventAction->fSiPMhitsX);
    analysisManager->CreateNtupleDColumn("SiPMHits_y_cm",
                                         fEventAction->fSiPMhitsY);
    analysisManager->CreateNtupleDColumn("SiPMHits_z_cm",
                                         fEventAction->fSiPMhitsZ);
    analysisManager->CreateNtupleDColumn("SiPMHits_Edep_keV",
                                         fEventAction->fSiPMhitsEdep);
    analysisManager->CreateNtupleDColumn(
        "SiPMHits_EdepNonIonizing_keV", fEventAction->fSiPMhitsEdepNonIonising);
    analysisManager->CreateNtupleDColumn("SiPMHits_TOA_ns",
                                         fEventAction->fSiPMhitsTOA);
    analysisManager->CreateNtupleIColumn("SiPMHits_type",
                                         fEventAction->fSiPMhitsType);

    analysisManager->CreateNtupleDColumn(
        "signalSum_HGCAL_GeV");                            // column Id = 23
    analysisManager->CreateNtupleDColumn("COGZ_HGCAL_cm"); // column Id = 24
    analysisManager->CreateNtupleIColumn("NHits_HGCAL");   // column Id = 25

    analysisManager->CreateNtupleDColumn(
        "signalSum_AHCAL_GeV");                            // column Id = 26
    analysisManager->CreateNtupleDColumn("COGZ_AHCAL_cm"); // column Id = 27
    analysisManager->CreateNtupleIColumn("NHits_AHCAL");   // column Id = 28
    analysisManager->FinishNtuple();
  }

  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *) {
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
