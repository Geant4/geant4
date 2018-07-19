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
/// \file electromagnetic/TestEm12/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 103638 2017-04-20 08:30:25Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "StepMax.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4EmCalculator.hh"
#include "G4EmParameters.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin)
:G4UserRunAction(),
 fDetector(det),fPhysics(phys),fPrimary(kin),fRun(0),fHistoManager(0)
{
  // Book predefined histograms
  fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fHistoManager; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector); 
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{    
  // show Rndm status
  if (isMaster) { 
    G4Random::showEngineStatus();
    G4EmParameters::Instance()->Dump();
  } 
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }
  
  if (!fPrimary) return;
      
  // keep run condition
  G4ParticleDefinition* particle 
      = fPrimary->GetParticleGun()->GetParticleDefinition();
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  fRun->SetPrimary(particle, energy);
     
  //get CsdaRange from EmCalculator
  //
  G4EmCalculator emCalculator;
  G4Material* material = fDetector->GetAbsorMaterial();

  G4double csdaRange = DBL_MAX;
  if (particle->GetPDGCharge() != 0.) {
    csdaRange = emCalculator.GetCSDARange(energy,particle,material);
    fRun->SetCsdaRange(csdaRange);
  }
     
  //set StepMax from histos 1 and 8
  //
  G4double stepMax = csdaRange;
  G4int ih = 1;
  if (analysisManager->GetH1Activation(ih))
    stepMax = analysisManager->GetH1Width(ih)*analysisManager->GetH1Unit(ih);
  ih = 8;
  if (analysisManager->GetH1Activation(ih)) {
    G4double width = analysisManager->GetH1Width(ih);
    stepMax = std::min(stepMax, width*csdaRange);        
  }                                        
  fPhysics->GetStepMaxProcess()->SetMaxStep2(stepMax);
  
  ///G4cout << "\n---> stepMax from histos 1 and 8 = " 
  ///       << G4BestUnit(stepMax,"Length") << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if (isMaster) fRun->EndOfRun(); 
  
  // save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {  
    analysisManager->Write();
    analysisManager->CloseFile();
  }      
 
  // show Rndm status
  if (isMaster) G4Random::showEngineStatus();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
