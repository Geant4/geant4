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
/// \file XrayTESdetHistoManager.cc
/// \brief Implementation of the HistoManager class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XrayTESdetHistoManager.hh"
#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayTESdetHistoManager::XrayTESdetHistoManager()
  : fFileName("example")
{
  fBook();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetHistoManager::fBook()
{
  // Create or get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetNtupleMerging(true);

  // Creating Ntuple
  int fNtupleID = analysisManager->CreateNtuple("TES_Tuple", "detectorPixelInfo");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleSColumn("vol_name");
  analysisManager->CreateNtupleIColumn("trackID");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("theta");
  analysisManager->CreateNtupleDColumn("phi");
  analysisManager->CreateNtupleIColumn("parentID");
  analysisManager->CreateNtupleIColumn("pixel_number");
  analysisManager->CreateNtupleDColumn("step_energy_dep");
  analysisManager->CreateNtupleIColumn("step_number");
  analysisManager->CreateNtupleDColumn("init_kinetic_energy");
  analysisManager->CreateNtupleDColumn("kinetic_energy");
  analysisManager->CreateNtupleSColumn("particle_name");
  analysisManager->CreateNtupleSColumn("pre_step_name");
  analysisManager->CreateNtupleSColumn("post_step_name");
  analysisManager->CreateNtupleSColumn("creator_proc_name");
  analysisManager->FinishNtuple(fNtupleID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
