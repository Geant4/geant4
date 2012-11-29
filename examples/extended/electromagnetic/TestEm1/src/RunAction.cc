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
/// \file electromagnetic/TestEm1/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:fDetector(det), fPrimary(kin)
{ 
 fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
 delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  ////G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  fNbOfTraks0 = fNbOfTraks1 = fNbOfSteps0 = fNbOfSteps1 = 0;
  fEdep = 0.;
  fTrueRange = fTrueRange2 = 0.;
  fProjRange = fProjRange2 = 0.;
  fTransvDev = fTransvDev2 = 0.;    
     
  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{ 
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  G4double dNbOfEvents = double(NbOfEvents);
    
  G4ParticleDefinition* particle = fPrimary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();    
  G4double energy   = fPrimary->GetParticleGun()->GetParticleEnergy();
  
  G4double length      = fDetector->GetSize();    
  G4Material* material = fDetector->GetMaterial();
  G4double density     = material->GetDensity();
   
  G4cout << "\n ======================== run summary ======================\n";
  
  G4int prec = G4cout.precision(5);
  
  G4cout << "\n The run was: " << NbOfEvents << " " << partName << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
         
 G4cout << "\n ============================================================\n";
 
 G4cout << "\n total energy deposit: " 
        << G4BestUnit(fEdep/dNbOfEvents, "Energy") << G4endl;
             
 //nb of tracks and steps per event
 //           
 G4cout << "\n nb tracks/event"
        << "   neutral: " << std::setw(10) << fNbOfTraks0/dNbOfEvents
        << "   charged: " << std::setw(10) << fNbOfTraks1/dNbOfEvents
        << "\n nb  steps/event"
        << "   neutral: " << std::setw(10) << fNbOfSteps0/dNbOfEvents
        << "   charged: " << std::setw(10) << fNbOfSteps1/dNbOfEvents
        << G4endl;
      
 //frequency of processes call
 std::map<G4String,G4int>::iterator it;         
 G4cout << "\n nb of process calls per event: \n   ";        
 for (it = fProcCounter.begin(); it != fProcCounter.end(); it++)  
    G4cout << std::setw(12) << it->first;
           
 G4cout << "\n   ";       
 for (it = fProcCounter.begin(); it != fProcCounter.end(); it++)   
    G4cout << std::setw(12) << (it->second)/dNbOfEvents;
 G4cout << G4endl;
      
 //compute true and projected ranges, and transverse dispersion
 //
 fTrueRange /= NbOfEvents; fTrueRange2 /= NbOfEvents;
 G4double trueRms = fTrueRange2 - fTrueRange*fTrueRange;        
 if (trueRms>0.) trueRms = std::sqrt(trueRms); else trueRms = 0.;
      
 fProjRange /= NbOfEvents; fProjRange2 /= NbOfEvents;
 G4double projRms = fProjRange2 - fProjRange*fProjRange;        
 if (projRms>0.) projRms = std::sqrt(projRms); else projRms = 0.;
       
 fTransvDev /= 2*NbOfEvents; fTransvDev2 /= 2*NbOfEvents;
 G4double trvsRms = fTransvDev2 - fTransvDev*fTransvDev;        
 if (trvsRms>0.) trvsRms = std::sqrt(trvsRms); else trvsRms = 0.;
 
 //compare true range with csda range from PhysicsTables
 //
 G4EmCalculator emCalculator;
 G4double rangeTable = 0.;
 if (particle->GetPDGCharge() != 0.)
   rangeTable = emCalculator.GetCSDARange(energy,particle,material);
      
 G4cout << "\n---------------------------------------------------------\n";
 G4cout << " Primary particle : " ;
 G4cout << "\n true Range = " << G4BestUnit(fTrueRange,"Length")
        << "   rms = "        << G4BestUnit(trueRms,  "Length");

 G4cout << "\n proj Range = " << G4BestUnit(fProjRange,"Length")
        << "   rms = "        << G4BestUnit(projRms,  "Length");
             
 G4cout << "\n proj/true  = " << fProjRange/fTrueRange;
                   
 G4cout << "\n transverse dispersion at end = " 
        << G4BestUnit(trvsRms,"Length");
        
 G4cout << "\n      mass true Range from simulation = " 
        << G4BestUnit(fTrueRange*density, "Mass/Surface")
        << "\n       from PhysicsTable (csda range) = " 
        << G4BestUnit(rangeTable*density, "Mass/Surface");        
 G4cout << "\n---------------------------------------------------------\n";
 G4cout << G4endl;
 
 // reset default precision
 G4cout.precision(prec);
                                    
  // remove all contents in fProcCounter 
  fProcCounter.clear();
    
  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }    
  
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
