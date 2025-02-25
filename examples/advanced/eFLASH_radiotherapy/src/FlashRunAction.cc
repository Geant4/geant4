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
/// \file FlashRunAction.cc
/// \brief Implementation of the FlashRunAction class

#include "FlashRunAction.hh"

#include "FlashPrimaryGeneratorAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4AnalysisManager.hh"

FlashRunAction::FlashRunAction() : G4UserRunAction()
{
  // Create analysis manager
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  
  // Merging ntuples from multiple root files due to multithreading 
  // fAnalysisManager->SetNtupleMerging(true);

  // Creation of ntuple
  fAnalysisManager->CreateNtuple("fSensitiveDetector", "fSensitiveDetector");
  fAnalysisManager->CreateNtupleDColumn("fX");
  fAnalysisManager->CreateNtupleDColumn("fY");
  fAnalysisManager->CreateNtupleDColumn("fZ");
  fAnalysisManager->CreateNtupleDColumn("fDose");
  fAnalysisManager->CreateNtupleDColumn("fEdep");
  fAnalysisManager->CreateNtupleIColumn("fEventID");
  fAnalysisManager->CreateNtupleIColumn("fParentID");
  fAnalysisManager->CreateNtupleSColumn("fParticle");
  fAnalysisManager->FinishNtuple(0);


}

FlashRunAction::~FlashRunAction() {}

void FlashRunAction::BeginOfRunAction(const G4Run *run) {
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  // Create analysis manager
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  
  // Merging ntuples from multiple root files due to multithreading
  //fAnalysisManager->SetNtupleMerging(true); // uncomment this only for root files

  // Open an output file
  G4int runID = run -> GetRunID();
  std::stringstream strRunID;
  strRunID << runID;
  fAnalysisManager->OpenFile("output_"+strRunID.str()+".csv");
  

  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void FlashRunAction::EndOfRunAction(const G4Run *run) {
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0)
    return;

  // close the output file
  G4AnalysisManager* fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();
  

  if (IsMaster()) {
    G4cout << G4endl
           << "--------------------End of Global Run-----------------------"
           << G4endl << "  The run was " << nofEvents << " events ";
  } else {
    G4cout << G4endl
           << "--------------------End of Local Run------------------------"
           << G4endl << "  The run was " << nofEvents << " events ";
  }
  
}
