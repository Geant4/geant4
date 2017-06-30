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
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $ID$
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4Run.hh"
#include "TrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "Analysis.hh"
#include "G4Threading.hh"
#include "CommandLineParser.hh"
//#include "NeuronLoadDataFile.hh"
#include "Run.hh"
#include "Randomize.hh"
#include <iomanip>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim) 
: G4UserRunAction(), //fpTrackingAction(0), fInitialized(0),
 fDebug(false),
fDetector(det),fPrimary(prim),fRun(0) 
{
   //CreateHistogram();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{  
  //delete G4AnalysisManager::Instance();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector); 
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{

RunInitManager::Instance()->Initialize();

  // keep run condition
  if ( fPrimary ) { 
    G4ParticleDefinition* particle 
      = fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }

/*
  G4cout << "##### Create analysis manager " << "  " << this << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFirstHistoId(1);
//  if(!analysisManager->IsActive()) {return; }

  G4cout << "Using " << analysisManager->GetType() <<
      " analysis manager" << G4endl;
   // Open an output file
  analysisManager->OpenFile("neuronG4");
  G4cout << "\n----> Histogram file is opened in neuronG4" <<
      "." << analysisManager->GetFileType() << G4endl;
// -----------------------------------------------------
  //Declare ntuples
  //
  // Create 1st ntuple (id = 1)
  //
  analysisManager->CreateNtuple("ntuple0", "Soma3D");
  analysisManager->CreateNtupleDColumn("ind");
  analysisManager->CreateNtupleDColumn("dist");
  analysisManager->CreateNtupleDColumn("edep");
  analysisManager->CreateNtupleDColumn("dose");
  analysisManager->FinishNtuple();
  //G4cout << "Ntuple-1 created" << G4endl;

  // Create 2nd ntuple (id = 2)
  //
  analysisManager->CreateNtuple("ntuple1", "Dend3D");
  analysisManager->CreateNtupleDColumn("indD");
  analysisManager->CreateNtupleDColumn("distD");
  analysisManager->CreateNtupleDColumn("edepD");
  analysisManager->CreateNtupleDColumn("doseD");
  analysisManager->FinishNtuple();
  //G4cout << "Ntuple-2 created" << G4endl;

  // Create 3rd ntuple (id = 3)
  //
  analysisManager->CreateNtuple("ntuple2", "Axon3D");
  analysisManager->CreateNtupleDColumn("indA");
  analysisManager->CreateNtupleDColumn("distA");
  analysisManager->CreateNtupleDColumn("edepA");
  analysisManager->CreateNtupleDColumn("doseA");
  analysisManager->FinishNtuple();
  //G4cout << "Ntuple-3 created" << G4endl;

  // Create 4rd ntuple (id = 4)
  //
  analysisManager->CreateNtuple("ntuple3", "Outputs per event");
  analysisManager->CreateNtupleDColumn("EdepAll");
  analysisManager->CreateNtupleDColumn("EdepMed");
  analysisManager->CreateNtupleDColumn("EdepSlice");
  analysisManager->CreateNtupleDColumn("EdepNeuron");
  analysisManager->FinishNtuple();
  //G4cout << "Ntuple-4 created" << G4endl;
  // ............................
  G4cout << "All Ntuples have been created " << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{ 

  if (isMaster) fRun->EndOfRun(); 

// save histogramms
/*  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
  //save histograms   
  analysisManager->Write();
  analysisManager->CloseFile();  
  // Complete clean-up
  delete G4AnalysisManager::Instance();
*/
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CreateHistogram()
{
  // Book histograms, ntuple

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in Analysis.hh

  CommandLineParser* parser = CommandLineParser::GetParser();
  Command* command(0);
  if((command = parser->GetCommandIfActive("-out"))==0) return;
//
// Declare ntuples
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::WriteHistogram()
{
  CommandLineParser* parser = CommandLineParser::GetParser();
  Command* commandLine(0);
  if((commandLine = parser->GetCommandIfActive("-out"))==0) return;

  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // if(!analysisManager->IsActive()) {return; }

  // save histograms
  //
  analysisManager->Write();
  analysisManager->CloseFile();

  if(fDebug)
  {
    G4cout << "================ ROOT FILES HAVE BEEN WRITTEN"
           << G4endl;
  }
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