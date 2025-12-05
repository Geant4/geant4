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
/// \file GB03RunAction.cc
/// \brief Implementation of the GB03RunAction class

#include "GB03RunAction.hh"

#include "GB03DetectorConstruction.hh"
#include "GB03Run.hh"

#include "G4AnalysisManager.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* GB03RunAction::GenerateRun()
{
  // Generate new RUN object
  return new GB03Run();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
GB03RunAction::GB03RunAction()
{
  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName("GB03");

  //
  // Book histograms
  //

  analysisManager->CreateH1("En", "E kin neutron", 100, 0., 100.0 * MeV);
  analysisManager->CreateH1("Eg", "E kin gamma", 100, 0., 15.0 * MeV);
  analysisManager->CreateH1("Nnhit", "The number of neutrons along x", 100, -1.5 * m, 150 * m / cm);
  analysisManager->CreateH1("Nghit", "The number of gammas along x", 100, -1.5 * m / cm,
                            1.5 * m / cm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  auto rm = G4RunManager::GetRunManager();
  auto det = static_cast<const GB03DetectorConstruction*>(rm->GetUserDetectorConstruction());
  fVerbose = det->GetVerboseLevel();

  auto* smeas = static_cast<G4Box*>(G4SolidStore::GetInstance()->GetSolid("Hodo"));

  auto halfX = smeas->GetXHalfLength() / cm;

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // set histo binning
  analysisManager->SetH1(2, 100, -halfX, halfX);
  analysisManager->SetH1(3, 100, -halfX, halfX);

  // Open an output file

  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03RunAction::EndOfRunAction(const G4Run* brun)
{
  auto run = static_cast<const GB03Run*>(brun);
  if (isMaster) {
    if (fVerbose > 0) {
      G4cout << "for the entire run " << run->GetNumberOfEvent() << G4endl << G4endl;
      run->GetPartCounter()->Print();
    }
  }
  else {
    if (fVerbose > 2) {
      G4cout << "for the local thread " << G4endl << G4endl;
      run->GetPartCounter()->Print();
    }
  }

  // Save Histograms
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
