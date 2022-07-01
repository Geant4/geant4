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
#include "HistoManager.hh"
#include "TrackingMessenger.hh"

#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction()
:G4UserTrackingAction(), fTrackMessenger(nullptr),
 fParticleCount(true)
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
	  
  G4int iabs        = track->GetTouchableHandle()->GetCopyNumber();
  G4String name     = track->GetDefinition()->GetParticleName();
  G4double meanLife = track->GetDefinition()->GetPDGLifeTime();  
  G4double energy   = track->GetKineticEnergy();
  if (fParticleCount && (iabs > 0))
    run->ParticleCount(iabs,name,energy,meanLife);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
  // get Run
  Run* run 
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  // where are we ?
  G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
  
  //status of primary particle : absorbed, transmited, reflected ?
  if (track->GetTrackID() == 1) {
    G4int flag = 0;
    if (status == fWorldBoundary) {
      if (track->GetMomentumDirection().x() > 0.) flag = 1;
      else                                        flag = 2;
    }
    run->AddTrackStatus(flag);
  }
  
  // keep only emerging particles
  if (status != fWorldBoundary) return;

  // count particles
  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4String name     = particle->GetParticleName();
  G4double meanLife = particle->GetPDGLifeTime();  
  G4double energy   = track->GetKineticEnergy();
  run->ParticleCount(0,name,energy,meanLife);

 ////G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

