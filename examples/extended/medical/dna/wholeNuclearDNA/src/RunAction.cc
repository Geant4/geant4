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
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Analysis.hh"
#include "G4Threading.hh"
#include "CommandLineParser.hh"

using namespace G4DNAPARSER;

void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() : G4UserRunAction(),
      fInitialized(0), fDebug(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
  // In this example, we considered that the same class was
  // used for both master and worker threads.
  // However, in case the run action is long,
  // for better code review, this practice is not recommanded.
  //
  // Please note, in the example provided with the Geant4 X beta version,
  // this RunAction class were not used by the master thread.


  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                     G4RunManager::sequentialRM);

  if(isMaster && sequential == false )
  // WARNING : in sequential mode, isMaster == true
  {
    BeginMaster(run);
  }
  else BeginWorker(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{

  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                     G4RunManager::sequentialRM);

  if(isMaster && sequential == false)
  {
    EndMaster(run);
  }
  else
  {
    EndWorker(run);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginMaster(const G4Run* run)
{
  if(fDebug)
  {
    bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                       G4RunManager::sequentialRM);
    G4cout << "===================================" << G4endl;
    if(!sequential)
      G4cout << "================ RunAction::BeginMaster" << G4endl;
    PrintRunInfo(run);
    G4cout << "===================================" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginWorker(const G4Run* run)
{
  if (fDebug)
  {
    G4cout << "===================================" << G4endl;
    G4cout << "================ RunAction::BeginWorker" << G4endl;
    PrintRunInfo(run);
    G4cout << "===================================" << G4endl;
  }
  if(fInitialized == false) InitializeWorker(run);

  CreateNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndMaster(const G4Run*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndWorker(const G4Run* run)
{
  if(fDebug)
  {
    PrintRunInfo(run);
  }

  G4int nofEvents = run->GetNumberOfEvent();
  if ( nofEvents == 0 )
  {
    if(fDebug)
    {
      G4cout << "================ NO EVENTS TREATED IN THIS RUN ==> Exit"
             << G4endl;
    }
    return;
  }

  ///////////////
  // Write Ntuple
  //
  WriteNtuple();

  ///////////////
  // Complete cleanup
  //
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::InitializeWorker(const G4Run*)
{
  // Initialize worker here
  // example: you want to retrieve pointers to other user actions
  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CreateNtuple()
{
  // Book histograms, ntuple

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in Analysis.hh

  CommandLineParser* parser = CommandLineParser::GetParser();
  Command* command(0);
  if((command = parser->GetCommandIfActive("-out"))==0) return;

  G4cout << "##### Create analysis manager " << "  " << this << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//  if(!analysisManager->IsActive()) {return; }

  G4cout << "Using " << analysisManager->GetType() <<
      " analysis manager" << G4endl;

  // Create directories

  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);

  // Open an output file
  G4String fileName;
  if(command->GetOption().empty() == false)
  {
    fileName = command->GetOption();
  }
  else
  {
   fileName = "wholeNuclearDNA";
//   fileName = command->GetDefaultOption(); // should work as well
  }
  analysisManager->OpenFile(fileName);

  // Creating ntuple
  analysisManager->CreateNtuple("ntuple", "geom_dna");
  analysisManager->CreateNtupleDColumn("flagParticle");
  analysisManager->CreateNtupleDColumn("flagProcess");
  analysisManager->CreateNtupleDColumn("flagVolume");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("edep");
  analysisManager->CreateNtupleDColumn("stepLength");

  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::WriteNtuple()
{
  CommandLineParser* parser = CommandLineParser::GetParser();
  Command* commandLine(0);
  if((commandLine = parser->GetCommandIfActive("-out"))==0) return;

  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//  if(!analysisManager->IsActive()) {return; }

  // save histograms
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintRunInfo(const G4Run* run)
{
  G4cout << "================ Run is = "
         << run->GetRunID() << G4endl;
  G4cout << "================ Run type is = "
         << G4RunManager::GetRunManager()->GetRunManagerType() << G4endl;
  G4cout << "================ Event processed = "
         << run->GetNumberOfEventToBeProcessed() << G4endl;
  G4cout << "================ Nevent = "
         << run->GetNumberOfEvent() << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
