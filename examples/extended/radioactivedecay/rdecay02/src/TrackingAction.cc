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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4HadronicProcessType.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det)
:G4UserTrackingAction(), fDetector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //which volume ?
  //
  G4LogicalVolume* lVolume = track->GetVolume()->GetLogicalVolume();
  G4int iVol = 0;
  if (lVolume == fDetector->GetLogicTarget())   iVol = 1;
  if (lVolume == fDetector->GetLogicDetector()) iVol = 2;
    
  //secondary particles only
  if (track->GetTrackID() == 1) return;
  
  const G4ParticleDefinition* particle = track->GetParticleDefinition();  
  G4String name   = particle->GetParticleName();
  G4int pid       = particle->GetPDGEncoding();
  G4int Z         = particle->GetAtomicNumber();
  G4int A         = particle->GetAtomicMass();
  G4double charge = particle->GetPDGCharge();    
  G4double energy = track->GetKineticEnergy();
  G4double time   = track->GetGlobalTime();
  G4double weight = track->GetWeight();
  
  run->ParticleCount(name,energy,iVol);
  
  //Radioactive decay products
  G4int processType = track->GetCreatorProcess()->GetProcessSubType();
  if (processType == fRadioactiveDecay) {
    //fill ntuple id = 3
    G4int id = 3;
    analysisManager->FillNtupleDColumn(id,0, double(pid));
    analysisManager->FillNtupleDColumn(id,1, double(Z));
    analysisManager->FillNtupleDColumn(id,2, double(A));
    analysisManager->FillNtupleDColumn(id,3, energy);
    analysisManager->FillNtupleDColumn(id,4, time/s);
    analysisManager->FillNtupleDColumn(id,5, weight);
    analysisManager->AddNtupleRow(id);
    
    if (charge < 3.) {   
      //fill ntuple id = 0
      id = 0;
      analysisManager->FillNtupleDColumn(id,0, double(pid));
      analysisManager->FillNtupleDColumn(id,1, energy);
      analysisManager->FillNtupleDColumn(id,2, time/s);
      analysisManager->FillNtupleDColumn(id,3, weight);
      analysisManager->AddNtupleRow(id);
    
      analysisManager->FillH1(6, energy, weight);
      analysisManager->FillH1(7, energy, weight);
      analysisManager->FillH1(8, energy, weight);
    }                        
  }
  
  //all unstable ions produced in target
  G4bool unstableIon = ((charge > 2.) && !(particle->GetPDGStable()));
  if ((unstableIon) && (iVol == 1)) {
    //fill ntuple id = 1
    G4int id = 1;
    analysisManager->FillNtupleDColumn(id,0, double(pid));
    analysisManager->FillNtupleDColumn(id,1, time/s);
    analysisManager->FillNtupleDColumn(id,2, weight);
    analysisManager->AddNtupleRow(id);  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* )
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

