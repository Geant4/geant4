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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4RunManager.hh"
#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction()
{ 
    // set printing event number per each 10 events
    //G4RunManager::GetRunManager()->SetPrintProgress(10);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{ 
    // inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    // book histograms, ntuple

    // create analysis manager
    // the choice of analysis technology is done via selection of a namespace
    // in Analysis.hh

    G4cout << "##### Create analysis manager " << "  " << this << G4endl;
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    G4cout << "Using " << analysisManager->GetType()
           << " analysis manager"
           << G4endl;

    analysisManager->SetVerboseLevel(1);

    // open an output file
    G4String fileName = "ICSD";
    analysisManager->OpenFile(fileName);

    // creating ntuple
    analysisManager->SetFirstNtupleId(1);

    // Creating ntuple
    analysisManager->CreateNtuple("ntuple_1", "ICSD");
    analysisManager->CreateNtupleDColumn(1,"ionisations");
    analysisManager->FinishNtuple(1);

    // Creating ntuple
    analysisManager->CreateNtuple("ntuple_2", "physics");
    analysisManager->CreateNtupleDColumn(2,"EventID");
    analysisManager->CreateNtupleDColumn(2,"flagProcess");
    analysisManager->CreateNtupleDColumn(2,"x");
    analysisManager->CreateNtupleDColumn(2,"y");
    analysisManager->CreateNtupleDColumn(2,"z");
    analysisManager->CreateNtupleDColumn(2,"totalEnergyDeposit");
    analysisManager->FinishNtuple(2);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* )
{
    // print histogram statistics

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    // save histograms

    analysisManager->Write();
    analysisManager->CloseFile();

    // complete cleanup

    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
