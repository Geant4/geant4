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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // Create analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // Define binning for each sort of histogram
  //
  //  Kinetic energy at the entrance and exit histograms
  const G4double kinEmin = 0.;
  const G4double kinEmax = 300.;
  const G4double kinEbinWidth = 0.1;
  const G4int kinEbins =
    static_cast<G4int>((kinEmax - kinEmin) / kinEbinWidth);

  // Logarithmic binning for energy histograms
  const G4int minLog10E = -4.;
  const G4int maxLog10E = +4.;
  const G4int nBinsLog10E = (maxLog10E - minLog10E) * 50;

  G4double binsLog10E[nBinsLog10E + 1];
  G4double binWidthLog10E =
    static_cast<G4double>((maxLog10E - minLog10E)) / nBinsLog10E;

  for (G4int ii = 0; ii <= nBinsLog10E; ii++) {
    binsLog10E[ii] = std::pow(10., binWidthLog10E * ii + minLog10E);
  }

  std::vector<G4double> vecBinsLog10E(binsLog10E, binsLog10E + nBinsLog10E + 1);

  // Logarithmic binning for weighted numbers histograms

  const G4int minLog10W = 0;
  const G4int maxLog10W = +6;
  const G4int nBinsLog10W = (maxLog10W - minLog10W) * 50;

  G4double binsLog10W[nBinsLog10W + 1];
  G4double binWidthLog10W =
    static_cast<G4double>((maxLog10W - minLog10W)) / nBinsLog10W;

  for (G4int ii = 0; ii <= nBinsLog10W; ii++) {
    binsLog10W[ii] = std::pow(10., binWidthLog10W * ii + minLog10W);
  }

  std::vector<G4double> vecBinsLog10W(binsLog10W, binsLog10W + nBinsLog10W + 1);

  // Linear binning for counter histograms
  G4int minCount = 0;
  G4int maxCount = 20000;
  G4int nBinsCount = (maxCount - minCount) / 100;

  // Create histograms
  //

  // Logarithmic binning  Energy histograms
  analysisManager->CreateH1(
    "fe",
    "Energy imparted per event [keV] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "efe",
    "Weighted energy imparted per event [keV] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "e2fe",
    "Squared-weighted energy imparted per event [keV] (log binning)",
    vecBinsLog10E);

  // Logarithmic binning for lineal energy histograms
  analysisManager->CreateH1("fy", "Lineal energy [keV/um] (log binning)",
                            vecBinsLog10E);

  analysisManager->CreateH1(
    "yfy",
    "Dose-weighted lineal energy [keV/um] (log binning)", vecBinsLog10E);

  analysisManager->CreateH1(
    "y2fy",
    "Squared-weighted lineal energy [keV/um] (log binning)", vecBinsLog10E);

  // Logarithmic binning for specific energy histograms
  analysisManager->CreateH1(
    "fz",
    "Single-event specific energy [Gy] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "zfz",
    "Dose-weighted single-event specific energy [Gy] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "z2fz",
    "Squared-weighted single-event specific energy [Gy] (log binning)",
    vecBinsLog10E);

  // Counters histograms
  analysisManager->CreateH1(
    "Nsel", "Number of selectable hits per event",
    nBinsCount, minCount, maxCount);

  analysisManager->CreateH1(
    "Nsite", "Number of hits in site", nBinsCount,
    minCount, maxCount);

  analysisManager->CreateH1("Nint",
    "Number of selectable hits in site", nBinsCount,
    minCount, maxCount);

  // Kinetic energy at the entrance and exit histograms
  analysisManager->CreateH1("KinE_in", "Kinetic energy at the entrance [MeV]",
                            kinEbins, kinEmin, kinEmax);
  analysisManager->CreateH1("KinE_out", "Kinetic energy at the exit [MeV]",
                            kinEbins, kinEmin, kinEmax);

  // 2D histogram for Nsite vs weighted Edep
  analysisManager->CreateH2(
    "Nsite_vs_e",
    "Number of hits in site vs energy imparted [keV] (log-log)",
    vecBinsLog10E, vecBinsLog10W);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
  // inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  // The file extension will set the choice of the output format
  G4String fileName = "microtrack.root";
  // Other formats supported: .csv, .hdf5, .xml
  // G4String fileName = "microtrack.csv";
  // G4String fileName = "microtrack.hdf5";
  // G4String fileName = "microtrack.xml";
  analysisManager->OpenFile(fileName);
  G4cout << "Using " << analysisManager->GetType() << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{
  auto analysisManager = G4AnalysisManager::Instance();

  if (IsMaster() && analysisManager->GetH1(1)) {

    G4cout << G4endl;

    G4cout << "----> print histogram statistics for the entire run: "
           << G4endl << G4endl;

    // print histogram statistics
    //

    G4cout << "  Single-event energy imparted:\n"
           << "  ----------------------------"
           << G4endl << G4endl;

    G4cout << "    Frequency-mean: \\varepsilon_{1,F} = "
           << analysisManager->GetH1(0)->mean() << " keV "
           << " (rms = " << analysisManager->GetH1(0)->rms() << " keV)"
           << G4endl;

    G4cout << "    Dose-mean: \\varepsilon_{1,D} = "
           << analysisManager->GetH1(1)->mean() << " keV "
           << " (rms = " << analysisManager->GetH1(1)->rms() << " keV)"
           << G4endl;

    G4cout << "    Squared-weighted histogram: mean = "
           << analysisManager->GetH1(2)->mean() << " keV "
           << " (rms = " << analysisManager->GetH1(2)->rms() << " keV)"
           << G4endl;

    G4cout << G4endl
           << "  Lineal energy:\n"
           << "  -------------"
           << G4endl << G4endl;

    G4cout << "    Frequency-mean: y_F = "
           << analysisManager->GetH1(3)->mean() << " keV/um "
           << " (rms = " << analysisManager->GetH1(3)->rms() << " keV/um)"
           << G4endl;

    G4cout << "    Dose-mean: y_D = "
           << analysisManager->GetH1(4)->mean() << " keV/um "
           << " (rms = " << analysisManager->GetH1(4)->rms() << " keV/um)"
           << G4endl;

    G4cout << "    Squared-weighted histogram: mean = "
           << analysisManager->GetH1(5)->mean() << " keV/um "
           << " (rms = " << analysisManager->GetH1(5)->rms() << " keV/um)"
           << G4endl;

    G4cout << G4endl
           << "  Single-event specific energy:\n"
           << "  ----------------------------"
           << G4endl << G4endl;

    G4cout << "    Frequency-mean: z_{1,F} = "
           << analysisManager->GetH1(6)->mean() << " Gy "
           << " (rms = " << analysisManager->GetH1(6)->rms() << " Gy)"
           << G4endl;

    G4cout << "    Dose-mean: z_{1,D} = "
           << analysisManager->GetH1(7)->mean() << " Gy "
           << " (rms = " << analysisManager->GetH1(7)->rms() << " Gy)"
           << G4endl;

    G4cout << "    Squared-weighted histogram: mean = "
           << analysisManager->GetH1(8)->mean() << " Gy"
           << " (rms = " << analysisManager->GetH1(8)->rms() << " Gy)"
           << G4endl;

    G4cout << G4endl
           << "  Number of hits per event:\n"
           << "  ------------------------"
           << G4endl << G4endl;

    G4cout << "    Eligible for site random placement, N_{sel}: mean = "
           << analysisManager->GetH1(9)->mean()
           << " (rms = " << analysisManager->GetH1(9)->rms() << ")"
           << G4endl;

    G4cout << "    Within the site, N_{site}: mean = "
           << analysisManager->GetH1(10)->mean()
           << " (rms = " << analysisManager->GetH1(10)->rms() << ")"
           << G4endl;

    G4cout << "    Within the site and eligible for site random placement,"
           << " N_{int}: mean = " << analysisManager->GetH1(11)->mean()
           << " (rms = " << analysisManager->GetH1(11)->rms() << ")"
           << G4endl;

    G4cout << G4endl
           << "  Kinetic energy of the primary particle:\n"
           << "  --------------------------------------"
           << G4endl << G4endl;

    G4cout << "    At entrance of SDbox, T_{in}: mean = "
           << G4BestUnit(analysisManager->GetH1(12)->mean(), "Energy")
           << " (rms = "
           << G4BestUnit(analysisManager->GetH1(12)->rms(), "Energy") << ")"
           << G4endl;

    G4cout << "    At exit of SDbox, T_{out}: mean = "
           << G4BestUnit(analysisManager->GetH1(13)->mean(), "Energy")
           << " (rms = "
           << G4BestUnit(analysisManager->GetH1(13)->rms(), "Energy") << ")"
           << G4endl;

    G4cout << G4endl;
  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
