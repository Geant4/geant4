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
/// \file electromagnetic/TestEm8/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//---------------------------------------------------------------------------
//
// ClassName:   RunAction
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "Run.hh"
#include "TestParameters.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
  : G4UserRunAction(), fAnalysisManager(0), fRun(0)
{
  TestParameters::GetPointer();
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->SetFileName("testem8");
  fAnalysisManager->SetVerboseLevel(1);
  fAnalysisManager->SetActivation(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Book()
{
  // Always creating analysis manager
  TestParameters* param = TestParameters::GetPointer();
  G4int nBinsE = param->GetNumberBins();
  G4int nBinsCluster = param->GetNumberBinsCluster();
  G4int nMaxCluster = param->GetMaxCluster();
  G4double maxEnergy = param->GetMaxEnergy();
  G4double factorALICE = param->GetFactorALICE();

  //G4cout << "### maxenergy(keV)= " << maxEnergy/keV 
  //         << "  factorALICE= " <<  factorALICE << G4endl;

  // Creating an 1-dimensional histograms in the root directory of the tree
  fAnalysisManager->SetFirstHistoId(1);   
  fAnalysisManager->CreateH1("h1","Energy deposition in detector (keV)",
                                  nBinsE,0.0,maxEnergy/keV);
  fAnalysisManager->CreateH1("h2","Number of primary clusters",
                                  nBinsCluster,0.0,G4double(nMaxCluster));
  fAnalysisManager->CreateH1("h3","Energy deposition in detector (ADC)",
                                  nBinsE,0.0,maxEnergy*factorALICE);
  fAnalysisManager->OpenFile(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(); 
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int id = aRun->GetRunID();
  G4cout << "### Run " << id << " start analysis activation: " 
         << G4endl;

  fRun->BeginOfRun();

  //histograms
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // print Run summary
  G4cout << "RunAction: End of run actions are started " << isMaster 
         << "  Nevt=  " << fRun->GetNumberOfEvent() 
         << "  Edep= " << fRun->GetStat()->mean() << G4endl;
  if (isMaster) {
    fRun->EndOfRun(); 
  }
  // save histos and close analysis
  if (fAnalysisManager->IsActive()) { 
    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();
  }
  //delete fAnalysisManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
