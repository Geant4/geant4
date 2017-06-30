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
/// \file exoticphysics/dmparticle/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 92047 2015-08-14 07:23:37Z gcosmo $
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
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
  : G4UserRunAction(), fAnalysisManager(nullptr), fRun(nullptr)
{
  TestParameters::GetPointer();
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->SetFileName("dmp");
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
  G4double maxEnergy = param->GetMaxEnergy();

  // Creating an 1-dimensional histograms in the root directory of the tree
  fAnalysisManager->SetFirstHistoId(1);   
  fAnalysisManager->CreateH1("h1","Energy deposition in detector (GeV)",
                             nBinsE,0.0,maxEnergy/GeV);
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

#ifdef G4VIS_USE
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
