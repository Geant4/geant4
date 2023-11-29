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
#include "Par04RunAction.hh"
#include <G4GenericAnalysisManager.hh>   // for G4GenericAnalysisManager
#include <G4ThreeVector.hh>              // for G4ThreeVector
#include <G4Types.hh>                    // for G4int, G4double
#include <G4UserRunAction.hh>            // for G4UserRunAction
#include "G4AnalysisManager.hh"          // for G4AnalysisManager
#include "Par04DetectorConstruction.hh"  // for Par04DetectorConstruction
#include "Par04EventAction.hh"           // for Par04EventAction

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04RunAction::Par04RunAction(Par04DetectorConstruction* aDetector, Par04EventAction* aEventAction)
  : G4UserRunAction()
  , fDetector(aDetector)
  , fEventAction(aEventAction)
{
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

  // Default filename, can be overriden with /analysis/setFileName
  analysisManager->SetFileName("Par04Output");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04RunAction::~Par04RunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04RunAction::BeginOfRunAction(const G4Run*)
{
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Create directories
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(0);

  // Get detector dimensions
  G4int cellNumZ       = fDetector->GetMeshNbOfCells().z();
  G4int cellNumRho     = fDetector->GetMeshNbOfCells().x();
  G4int cellNumPhi     = fDetector->GetMeshNbOfCells().y();
  G4double cellSizeZ   = fDetector->GetMeshSizeOfCells().z();
  G4double cellSizeRho = fDetector->GetMeshSizeOfCells().x();
  G4double cellSizePhi = fDetector->GetMeshSizeOfCells().y();
  // Default max value of energy stored in histogram (in GeV)
  G4double maxEnergy = 1000;

  // Creating control histograms
  analysisManager->CreateH1("energyParticle", "Primary energy;E_{MC} (GeV);Entries", 1024, 0,
                            1.1 * maxEnergy);
  analysisManager->CreateH1("energyDeposited", "Deposited energy;E_{MC} (GeV);Entries", 1024, 0,
                            1.1 * maxEnergy);
  analysisManager->CreateH1(
    "energyRatio", "Ratio of energy deposited to primary;E_{dep} /  E_{MC};Entries", 1024, 0, 1);
  analysisManager->CreateH1("time", "Simulation time; time (s);Entries", 2048, 0, 100);
  analysisManager->CreateH1("longProfile", "Longitudinal profile;t (mm);#LTE#GT (MeV)", cellNumZ,
                            -0.5 * cellSizeZ, (cellNumZ - 0.5) * cellSizeZ);
  analysisManager->CreateH1("transProfile", "Transverse profile;r (mm);#LTE#GT (MeV)", cellNumRho,
                            -0.5 * cellSizeRho, (cellNumRho - 0.5) * cellSizeRho);
  analysisManager->CreateH1("longFirstMoment",
                            "First moment of longitudinal distribution;#LT#lambda#GT (mm);Entries",
                            1024, -0.5 * cellSizeZ,
                            cellNumZ * cellSizeZ / 2);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1("transFirstMoment",
                            "First moment of transverse distribution;#LTr#GT "
                            "(mm);Entries",
                            1024, -0.5 * cellSizeRho,
                            cellNumRho * cellSizeRho /
                              1);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "longSecondMoment",
    "Second moment of longitudinal distribution;#LT#lambda^{2}#GT "
    "(mm^{2});Entries",
    1024, 0, std::pow(cellNumZ * cellSizeZ, 2) / 25);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1(
    "transSecondMoment", "Second moment of transverse distribution;#LTr^{2}#GT (mm^{2});Entries",
    1024, 0, std::pow(cellNumRho * cellSizeRho, 2) / 5);  // arbitrary scaling of max value on axis
  analysisManager->CreateH1("hitType", "hit type;type (0=full, 1= fast);Entries", 2, -0.5, 1.5);
  analysisManager->CreateH1("phiProfile", "Azimuthal angle profile, centred at mean;phi;#LTE#GT (MeV)", cellNumPhi,
                            - (cellNumPhi - 0.5) * cellSizePhi, (cellNumPhi - 0.5) * cellSizePhi);
  analysisManager->CreateH1("numHits", "Number of hits above 0.5 keV", 4048, 0, 40500);
  analysisManager->CreateH1("cellEnergy", "Cell energy distribution;log10(E/MeV);Entries", 1024, -4, 2);

  // Creating ntuple
  analysisManager->CreateNtuple("events", "per event data");
  analysisManager->CreateNtupleDColumn("EnergyMC");
  analysisManager->CreateNtupleDColumn("EnergyCell", fEventAction->GetCalEdep());
  analysisManager->CreateNtupleIColumn("rhoCell", fEventAction->GetCalRho());
  analysisManager->CreateNtupleIColumn("phiCell", fEventAction->GetCalPhi());
  analysisManager->CreateNtupleIColumn("zCell", fEventAction->GetCalZ());
  analysisManager->CreateNtupleDColumn("SimTime");
  analysisManager->FinishNtuple();

  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04RunAction::EndOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}
