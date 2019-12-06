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
/// \file B5RunAction.cc
/// \brief Implementation of the B5RunAction class

#include "B5RunAction.hh"
#include "B5EventAction.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "g4analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5RunAction::B5RunAction(B5EventAction* eventAction)
 : G4UserRunAction(),
   fEventAction(eventAction)
{ 
  // Create the analysis manager using a new factory method.
  // The choice of analysis technology is done via the function argument.
  auto analysisManager = G4Analysis::ManagerInstance("root");
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetNtupleMerging(true);
     // Note: merging ntuples is available only with Root output
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("B5");

  // Book histograms, ntuple
  //
  
  // Creating 1D histograms
  analysisManager
    ->CreateH1("Chamber1","Drift Chamber 1 # Hits", 50, 0., 50); // h1 Id = 0
  analysisManager
    ->CreateH1("Chamber2","Drift Chamber 2 # Hits", 50, 0., 50); // h1 Id = 1
  
  // Creating 2D histograms
  analysisManager                                                
    ->CreateH2("Chamber1 XY","Drift Chamber 1 X vs Y",           // h2 Id = 0
               50, -1000., 1000, 50, -300., 300.); 
  analysisManager                                                
    ->CreateH2("Chamber2 XY","Drift Chamber 2 X vs Y",           // h2 Id = 1
               50, -1500., 1500, 50, -300., 300.);

  // Creating ntuple
  //
  if ( fEventAction ) {
    analysisManager->CreateNtuple("B5", "Hits");
    analysisManager->CreateNtupleIColumn("Dc1Hits");  // column Id = 0
    analysisManager->CreateNtupleIColumn("Dc2Hits");  // column Id = 1
    analysisManager->CreateNtupleDColumn("ECEnergy"); // column Id = 2
    analysisManager->CreateNtupleDColumn("HCEnergy"); // column Id = 3
    analysisManager->CreateNtupleDColumn("Time1");    // column Id = 4
    analysisManager->CreateNtupleDColumn("Time2");    // column Id = 5
    analysisManager                                   // column Id = 6
      ->CreateNtupleDColumn("ECEnergyVector", fEventAction->GetEmCalEdep()); 
    analysisManager                                   // column Id = 7
      ->CreateNtupleDColumn("HCEnergyVector", fEventAction->GetHadCalEdep());
    analysisManager->FinishNtuple();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5RunAction::~B5RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file 
  // The default file name is set in B5RunAction::B5RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
