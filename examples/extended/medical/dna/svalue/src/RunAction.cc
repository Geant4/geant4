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
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "Run.hh"
#include "HistoManager.hh"
#include "MyFile.hh"

#ifdef MYFILE
 #include "MyPrimaryGeneratorActionFromFile.hh"
#else
 #include "PrimaryGeneratorAction.hh"
#endif

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

#ifdef MYFILE
   
   const MyPrimaryGeneratorActionFromFile* primary =
      dynamic_cast<const MyPrimaryGeneratorActionFromFile*>(G4RunManager::GetRunManager()
          ->GetUserPrimaryGeneratorAction());
   
#else
   const PrimaryGeneratorAction* primary =
      dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()
          ->GetUserPrimaryGeneratorAction());

#endif

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
