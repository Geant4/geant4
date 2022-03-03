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
#include "Par03RunAction.hh"
#include "Par03DetectorConstruction.hh"

#include "G4AnalysisManager.hh"

Par03RunAction::Par03RunAction(Par03DetectorConstruction* aDetector)
  : G4UserRunAction()
  , fDetector(aDetector)
{
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  // Default filename, can be overriden with /analysis/setFileName
  analysisManager->SetFileName("Par03Output");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03RunAction::~Par03RunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03RunAction::BeginOfRunAction(const G4Run*)
{
  // Get detector dimensions
  G4int cellNumZ       = fDetector->GetNbOfLayers();
  G4int cellNumRho     = fDetector->GetNbOfRhoCells();
  G4double cellSizeZ   = fDetector->GetLength() / cellNumZ;
  G4double cellSizeRho = fDetector->GetRadius() / cellNumRho;
  // Default max value of energy stored in histogram (in GeV)
  G4double maxEnergy = 100;

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Creating control histograms
  analysisManager->CreateH1("energyParticle",
                            "Primary energy;E_{MC} (GeV);Entries", 256, 0,
                            1.1 * maxEnergy);
  analysisManager->CreateH1("energyDeposited",
                            "Deposited energy;E_{MC} (GeV);Entries", 256, 0,
                            1.1 * maxEnergy);
  analysisManager->CreateH1(
    "energyRatio",
    "Ratio of energy deposited to primary;E_{dep} /  E_{MC};Entries", 1024, 0,
    1);
  analysisManager->CreateH1("time", "Simulation time; time (s);Entries", 2048,
                            0, 30);
  analysisManager->CreateH1(
    "longProfile", "Longitudinal profile;t (mm);#LTE#GT (MeV)", cellNumZ,
    -0.5 * cellSizeZ, (cellNumZ - 0.5) * cellSizeZ);
  analysisManager->CreateH1(
    "transProfile", "Transverse profile;r (mm);#LTE#GT (MeV)", cellNumRho,
    -0.5 * cellSizeRho, (cellNumRho - 0.5) * cellSizeRho);
  analysisManager->CreateH1(
    "longFirstMoment",
    "First moment of longitudinal distribution;#LT#lambda#GT (mm);Entries",
    1024, -0.5 * cellSizeZ,
    cellNumZ * cellSizeZ / 2);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1("transFirstMoment",
                            "First moment of transverse distribution;#LTr#GT "
                            "(mm);Entries",
                            1024, -0.5 * cellSizeRho,
                            cellNumRho * cellSizeRho /
                              10);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "longSecondMoment",
    "Second moment of longitudinal distribution;#LT#lambda^{2}#GT "
    "(mm^{2});Entries",
    1024, 0,
    std::pow(cellNumZ * cellSizeZ, 2) /
      25);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "transSecondMoment",
    "Second moment of transverse distribution;#LTr^{2}#GT (mm^{2});Entries",
    1024, 0,
    std::pow(cellNumRho * cellSizeRho, 2) /
      25);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "hitType", "hit type;type (0=full, 1= fast);Entries", 2, -0.5, 1.5);

  // Open an output file
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03RunAction::EndOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  analysisManager->Clear();
}