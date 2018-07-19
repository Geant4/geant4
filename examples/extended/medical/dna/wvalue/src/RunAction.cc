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
/// \file medical/dna/wvalue/src/RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
:G4UserRunAction(),
fpDetector(0), fpRun(0),fpHistoManager(0)
{
  fpDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction());

  // Book predefined histograms
  fpHistoManager = new HistoManager();
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fpHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fpRun = new Run(fpDetector);
  return fpRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{    
  // save Rndm status
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  // if (isMaster) G4Random::showEngineStatus();
  
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }
  
  //============================================================================
  const PrimaryGeneratorAction* primary =
      dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()
          ->GetUserPrimaryGeneratorAction());

  if (!primary) return; //
      
  // keep run condition
  G4ParticleDefinition* particle 
      = primary->GetParticleGun()->GetParticleDefinition();
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();
  fpRun->SetPrimary(particle, energy);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if (isMaster) fpRun->EndOfRun();
  
  // save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {  
    analysisManager->Write();
    analysisManager->CloseFile();
  }      
 
  // show Rndm status
  // if (isMaster) G4Random::showEngineStatus();  
}
