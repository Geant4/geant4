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

#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "TrackingMessenger.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4ParticleTypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(EventAction* event)
:G4UserTrackingAction(), fEventAction(event), fTrackMessenger(nullptr),
 fParticleCount(true) , fKillNeutron(false)
{   
  fTrackMessenger = new TrackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
  delete fTrackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{  
  //count secondary particles
  if (track->GetTrackID() == 1) return;
    
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
	
  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4String name     = particle->GetParticleName();
  G4double meanLife = particle->GetPDGLifeTime();
  G4double energy   = track->GetKineticEnergy();
  if (fParticleCount) run->ParticleCount(name,energy,meanLife);
       
  // histograms: energy at creation
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
 
  G4int ih = 0;
  G4String type   = particle->GetParticleType();      
  G4double charge = particle->GetPDGCharge();
  if (charge > 3.)  ih = 10; 
  else if (particle == G4Gamma::Gamma())       ih = 4;
  else if (particle == G4Electron::Electron()) ih = 5;
  else if (particle == G4Positron::Positron()) ih = 5;  
  else if (particle == G4Neutron::Neutron())   ih = 6;
  else if (particle == G4Proton::Proton())     ih = 7;
  else if (particle == G4Deuteron::Deuteron()) ih = 8;
  else if (particle == G4Alpha::Alpha())       ih = 9;       
  else if (type == "nucleus")                  ih = 10;
  else if (type == "baryon")                   ih = 11;         
  else if (type == "meson")                    ih = 12;
  else if (type == "lepton")                   ih = 13;        
  if (ih > 0) analysis->FillH1(ih,energy);
   
  //to force only 1 fission : kill secondary neutrons
  if (fKillNeutron && (particle == G4Neutron::Neutron())) {
    fEventAction->AddEdep(energy);  
    G4Track* aTrack = (G4Track*)track;
    aTrack->SetTrackStatus(fStopAndKill);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
 // keep only emerging particles
 G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
 if (status != fWorldBoundary) return; 
 
 const G4ParticleDefinition* particle = track->GetParticleDefinition();
 G4String name   = particle->GetParticleName();
 G4double energy = track->GetKineticEnergy();
 
 fEventAction->AddEflow(energy);  
 
 Run* run = static_cast<Run*>(
              G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 run->ParticleFlux(name,energy);               
 
 // histograms: energy at exit
 //
 G4AnalysisManager* analysis = G4AnalysisManager::Instance();
 
 G4int ih = 0; 
 G4String type   = particle->GetParticleType();      
 G4double charge = particle->GetPDGCharge();
 if (charge > 3.)  ih = 20; 
 else if (particle == G4Gamma::Gamma())       ih = 14;
 else if (particle == G4Electron::Electron()) ih = 15;
 else if (particle == G4Positron::Positron()) ih = 15;  
 else if (particle == G4Neutron::Neutron())   ih = 16;
 else if (particle == G4Proton::Proton())     ih = 17;
 else if (particle == G4Deuteron::Deuteron()) ih = 18;
 else if (particle == G4Alpha::Alpha())       ih = 19;       
 else if (type == "nucleus")                  ih = 20;
 else if (type == "baryon")                   ih = 21;         
 else if (type == "meson")                    ih = 22;
 else if (type == "lepton")                   ih = 23;        
 if (ih > 0) analysis->FillH1(ih,energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

