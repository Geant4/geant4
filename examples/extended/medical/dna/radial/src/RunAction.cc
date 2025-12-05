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

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det)
:fMyDetectorConstruction(det)
{
  if (isMaster)
  {
    G4cout << "##### Create analysis manager " << "  " << this << G4endl;
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->SetDefaultFileType("root");
    analysisManager->SetFirstNtupleId(1);

    // Create ntuple
    analysisManager->CreateNtuple("radial", "radial");
    analysisManager->CreateNtupleDColumn("radius");
    analysisManager->CreateNtupleDColumn("dose");
    analysisManager->FinishNtuple();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  fRun = new Run(fMyDetectorConstruction);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  if (isMaster)
  {
    auto analysisManager = G4AnalysisManager::Instance();

    // Open an output file
    G4String fileName = "radial";
    analysisManager->OpenFile(fileName);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if (isMaster)
  {
    // Display results from merged local runs
    fRun->EndOfRun();

    // Fill ntuple
    auto analysisManager = G4AnalysisManager::Instance();

    G4double cumulatedDeposit = 0;

    // Loop on cylinders and collect dose from merged local runs
    G4int nbCyl = fMyDetectorConstruction->GetCylinderNumber();

    for (G4int i = 0; i < nbCyl ; i++)
    {
      cumulatedDeposit = fRun->GetCylDoseDeposit(i);
      if (cumulatedDeposit > 0.)
      {
        analysisManager->FillNtupleDColumn
          (1,0,i*fMyDetectorConstruction->GetCylinderThickness()/nm);
        analysisManager->FillNtupleDColumn
          (1,1,cumulatedDeposit/fRun->GetNumberOfEvent()/gray);
        analysisManager->AddNtupleRow(1);
      }
    }

    // Save histograms
    analysisManager->Write();
    analysisManager->CloseFile();
    analysisManager->Clear();
  }
}
