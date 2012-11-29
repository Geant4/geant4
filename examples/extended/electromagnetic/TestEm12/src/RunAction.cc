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
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "StepMax.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin, HistoManager* histo)
:fDetector(det),fPhysics(phys),fKinematic(kin),fHistoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  ////G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
  
  //initialize total energy deposit
  //
  fEdeposit = fEdeposit2 = 0.; 
   
  //initialize track legth of primary
  //
  fTrackLen = fTrackLen2 = 0.;
    
  //initialize projected range
  //
  fProjRange = fProjRange2 = 0.;
    
  //initialize mean step size
  //
  fNbOfSteps = fNbOfSteps2 = 0;  fStepSize = fStepSize2 = 0.;
  
  //get fCsdaRange from EmCalculator
  //
  G4EmCalculator emCalculator;
  G4Material* material = fDetector->GetAbsorMaterial();
  G4ParticleDefinition* particle = fKinematic->GetParticleGun()
                                          ->GetParticleDefinition();
  G4double energy = fKinematic->GetParticleGun()->GetParticleEnergy();
  fCsdaRange = DBL_MAX;
  if (particle->GetPDGCharge() != 0.)
    fCsdaRange = emCalculator.GetCSDARange(energy,particle,material);
                         
  //histograms
  //
  fHistoManager->book();
  
  //set StepMax from histos
  //
  G4double stepMax = fHistoManager->ComputeStepMax(fCsdaRange);
  fPhysics->GetStepMaxProcess()->SetMaxStep2(stepMax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  
  G4int NbofEvents = aRun->GetNumberOfEvent();
  if (NbofEvents == 0) return;
  G4double fNbofEvents = double(NbofEvents);

  //run conditions
  //  
  G4Material* material = fDetector->GetAbsorMaterial();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = fKinematic->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();                             
  G4double energy = fKinematic->GetParticleGun()->GetParticleEnergy();
  
  G4cout << "\n ======================== run summary ======================\n";
  
  G4int prec = G4cout.precision(2);
  
  G4cout << "\n The run consists of " << NbofEvents << " "<< partName << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(fDetector->GetAbsorRadius(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
         
  G4cout << "\n ============================================================\n";
  
  //compute total energy deposit
  //
  fEdeposit /= NbofEvents; fEdeposit2 /= NbofEvents;
  G4double rms = fEdeposit2 - fEdeposit*fEdeposit;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Total Energy deposited        = " << G4BestUnit(fEdeposit,"Energy")
    << " +- "                               << G4BestUnit( rms,"Energy")
    << G4endl;
              
  //compute track length of primary track
  //
  fTrackLen /= NbofEvents; fTrackLen2 /= NbofEvents;
  rms = fTrackLen2 - fTrackLen*fTrackLen;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\n Track length of primary track = " << G4BestUnit(fTrackLen,"Length")
    << " +- "                               << G4BestUnit( rms,"Length");
    
  //compare with csda range
  //
  //G4EmCalculator emCalculator;
  //G4double fCsdaRange = 0.;
  //if (particle->GetPDGCharge() != 0.)
  //  fCsdaRange = emCalculator.GetCSDARange(energy,particle,material);
  G4cout 
    << "\n Range from EmCalculator       = " << G4BestUnit(fCsdaRange,"Length")
    << " (from full dE/dx)" << G4endl;
                     
  //compute projected range of primary track
  //
  fProjRange /= NbofEvents; fProjRange2 /= NbofEvents;
  rms = fProjRange2 - fProjRange*fProjRange;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;
   
  G4cout 
    << "\n Projected range               = " << G4BestUnit(fProjRange,"Length")
    << " +- "                                << G4BestUnit( rms,"Length")    
    << G4endl;
    
  //nb of steps and step size of primary track
  //
  G4double fNbSteps = fNbOfSteps/fNbofEvents, fNbSteps2 = fNbOfSteps2/fNbofEvents;
  rms = fNbSteps2 - fNbSteps*fNbSteps;       
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(2);       
  G4cout << "\n Nb of steps of primary track  = " << fNbSteps << " +- " << rms;
    
  fStepSize /= NbofEvents; fStepSize2 /= NbofEvents;
  rms = fStepSize2 - fStepSize*fStepSize;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(3);       
  G4cout 
    << "\t Step size= " << G4BestUnit(fStepSize,"Length")
    << " +- "           << G4BestUnit( rms,"Length")
    << G4endl;
     
  // normalize histogram of longitudinal energy profile
  //
  G4int ih = 1;
  G4double binWidth = fHistoManager->GetBinWidth(ih);
  G4double fac = (1./(NbofEvents*binWidth))*(mm/MeV);
  fHistoManager->Normalize(ih,fac);
  
  // normalize histogram d(E/E0)/d(r/r0)
  //
  ih = 8;
  binWidth = fHistoManager->GetBinWidth(ih);
  fac = 1./(NbofEvents*binWidth*energy);
  fHistoManager->Normalize(ih,fac);
    
   // reset default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  // save histograms
  fHistoManager->save();
 
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
