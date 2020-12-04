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
// Authors: S. Meylan and C. Villagrasa (IRSN, France)
// Updated: H. Tran (IRSN), France: 20/12/2018

#include "RunAction.hh"
#include "G4Run.hh"
#include "g4root.hh"
#include "globals.hh"
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction()
    : G4UserRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run*)
{
    CreateNtuple("output.root");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run*)
{
    WriteNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::CreateNtuple(const G4String& name)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetVerboseLevel(0);
    analysisManager->SetNtupleDirectoryName("ntuple");

    // open output file
    //
    G4bool fileOpen = analysisManager->OpenFile(name);
    if (!fileOpen)
    {
        G4cout << "\n---> HistoManager::book(): cannot open " << name<< G4endl;
        return;
    }

    analysisManager->SetFirstNtupleId(1);
    analysisManager->CreateNtuple("ntuple_1","physical_stage");
    analysisManager->CreateNtupleDColumn(1,"x");
    analysisManager->CreateNtupleDColumn(1,"y");
    analysisManager->CreateNtupleDColumn(1,"z");
    analysisManager->CreateNtupleDColumn(1,"edep");
    analysisManager->CreateNtupleDColumn(1,"diffKin");
    analysisManager->CreateNtupleIColumn(1,"volumeName");
    analysisManager->CreateNtupleDColumn(1,"CopyNumber");
    analysisManager->CreateNtupleIColumn(1,"EventID");
    analysisManager->FinishNtuple(1);

    // For chemistry

    analysisManager->CreateNtuple("ntuple_2","chem_stage");
    analysisManager->CreateNtupleDColumn(2,"x");
    analysisManager->CreateNtupleDColumn(2,"y");
    analysisManager->CreateNtupleDColumn(2,"z");
    analysisManager->CreateNtupleSColumn(2,"RadName");
    analysisManager->CreateNtupleIColumn(2,"EventID");
    analysisManager->FinishNtuple(2);

    G4cout << "\n----> Histogram file is opened in " << name << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void RunAction::WriteNtuple()
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    G4cout << "\n----> Histograms are saved" << G4endl;
}
