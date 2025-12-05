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
#include "G4Run.hh"
#include <CLHEP/Units/SystemOfUnits.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction()
{
    //using analysis manager for output
    auto analysisManager = G4AnalysisManager::Instance();

    //setting our histogram
    //a true range and bin number is set up in BeginOfRunAction
    analysisManager->CreateH1("Spectrum_electrons","Spectrum_e-",20,0,100);
    analysisManager->CreateH1("Spectrum_positrons","Spectrum_e+",20,0,100);
    analysisManager->CreateH1("Spectrum_gamma","Spectrum_gamma",25,0,125);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
    //opening output file
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4String fileName = "results.root";
    analysisManager->OpenFile(fileName);

    //setting histograms
    analysisManager->SetH1(0,20,0,100);
    analysisManager->SetH1XAxisTitle(0,"E_e- [GeV]");
    analysisManager->SetH1YAxisTitle(0,"Count");

    analysisManager->SetH1(1,20,0,100);
    analysisManager->SetH1XAxisTitle(1,"E_e+ [GeV]");
    analysisManager->SetH1YAxisTitle(1,"Count");

    analysisManager->SetH1(2,25,0,125);
    analysisManager->SetH1XAxisTitle(2,"E_gamma [GeV]");
    analysisManager->SetH1YAxisTitle(2,"Count");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
