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
/// \file electromagnetic/TestEm4/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);  
  analysisManager->SetFirstHistoId(1);
    
  // Creating histograms
  //
  analysisManager->CreateH1("1","energy (MeV) deposited in C6F6",100,0.,10.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();

  
   // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "testem4";
  analysisManager->OpenFile(fileName);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (isMaster) {
    // Run conditions
    //  note: There is no primary generator action object for "master"
    //        run manager for multi-threaded mode.
    PrimaryGeneratorAction* generatorAction = (PrimaryGeneratorAction*)
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    G4String runCondition = " ";
    if (generatorAction)
    {
      G4ParticleGun* particleGun = generatorAction->GetParticleGun();
      runCondition  = particleGun->GetParticleDefinition()->GetParticleName();
      runCondition += " of ";
      G4double particleEnergy = particleGun->GetParticleEnergy();
      runCondition += G4BestUnit(particleEnergy,"Energy");
      runCondition += " in a volume of C6F6 ";
    }

    G4cout
     << G4endl
     << "--------------------End of Global Run------------------------";
    G4cout
     << G4endl
     << " The run consists of " << nofEvents << " " << runCondition << G4endl;
    if ( analysisManager->GetH1(1) ) {
      G4cout << " Energy deposit per event:  mean = "
             << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
             << "    rms = "
             << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy")
	     << G4endl
             << "------------------------------------------------------------"
             << G4endl;
    }
  }

  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();

  //save histograms      
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
