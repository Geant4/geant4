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
/// \file electromagnetic/TestEm11/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 76258 2013-11-08 11:36:51Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "PhysicsList.hh"
#include "StepMax.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det),fPhysics(phys),fKinematic(kin),fRun(0),
 fHistoManager(0)
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
  // save Rndm status
  ////G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4Random::showEngineStatus();
  
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  } 
    
  //set StepMax from histos 1 and 8
  //e
  G4double stepMax = DBL_MAX;
  G4int ih = 1;
  if (analysisManager->GetH1Activation(ih)) {
    stepMax = analysisManager->GetH1Width(ih);
  }  
  
  if (  ! fKinematic ) return;
  
    //
  ih = 8;
  G4ParticleDefinition* particle = fKinematic->GetParticleGun()
                                          ->GetParticleDefinition();
  if (particle->GetPDGCharge() != 0.) {
    G4double width = analysisManager->GetH1Width(ih);
    if (width == 0.) width = 1.;                                               
    G4EmCalculator emCalculator;
    G4double energy = fKinematic->GetParticleGun()->GetParticleEnergy();
    G4int nbOfAbsor = fDetector->GetNbOfAbsor();
    for (G4int i=1; i<= nbOfAbsor; i++) {
      G4Material* material = fDetector->GetAbsorMaterial(i);
      G4double newCsdaRange 
                   = emCalculator.GetCSDARange(energy,particle,material);
      fRun->SetCsdaRange(i, newCsdaRange);
      if (analysisManager->GetH1Activation(ih))
        stepMax = std::min(stepMax, width*newCsdaRange);
      if (i>1) {
        G4double thickness = fDetector->GetAbsorThickness(i-1);
        G4double xfrontNorm = fRun->GetXfrontNorm(i-1);
        G4double csdaRange = fRun->GetCsdaRange(i-1);
        G4double newXfrontNorm = xfrontNorm + thickness/csdaRange;
        fRun->SetXfrontNorm(i, newXfrontNorm);
      }                              
    }        
  }    
  fPhysics->GetStepMaxProcess()->SetMaxStep2(stepMax);           
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(2);
 
  G4int nbofEvents = run->GetNumberOfEvent();

  if ( fKinematic && nbofEvents ) {

    //run conditions
    //     
    G4ParticleDefinition* particle = fKinematic->GetParticleGun()
                                          ->GetParticleDefinition();
    G4String partName = particle->GetParticleName();
    G4double energy = fKinematic->GetParticleGun()->GetParticleEnergy();
  
    G4int nbOfAbsor = fDetector->GetNbOfAbsor();
  
    G4cout << "\n ======================== run summary =====================\n";
  
    G4cout 
      << "\n The run consists of " << nbofEvents << " "<< partName << " of "
      << G4BestUnit(energy,"Energy") 
      << " through "  << nbOfAbsor << " absorbers: \n";
    for (G4int i=1; i<= nbOfAbsor; i++) {
       G4Material* material = fDetector->GetAbsorMaterial(i);
       G4double thickness = fDetector->GetAbsorThickness(i);
       G4double density = material->GetDensity();    
       G4cout << std::setw(20) << G4BestUnit(thickness,"Length") << " of "
              << material->GetName() << " (density: " 
              << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
    }         
    G4cout << "\n ==========================================================\n";
  }
  
  //compute and print total energy deposit
  //
  if ( isMaster ) {
    fRun->ComputeStatistics();  

    // normalize histograms of longitudinal energy profile
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int ih = 1;
    G4double binWidth = analysisManager->GetH1Width(ih);
    G4double fac = (1./(nbofEvents*binWidth))*(mm/MeV);
    analysisManager->ScaleH1(ih,fac);
    
    ih = 8;
    binWidth = analysisManager->GetH1Width(ih);
    fac = (1./(nbofEvents*binWidth))*(g/(MeV*cm2));
    analysisManager->ScaleH1(ih,fac);
  }

  // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  // save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {  
    analysisManager->Write();
    analysisManager->CloseFile();
  }    
 
  // show Rndm status
  G4Random::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
