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
#include "G4ios.hh"

#include "exrdmSteppingAction.hh"
#include "exrdmAnalysisManager.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmSteppingAction::exrdmSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdmSteppingAction::~exrdmSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdmSteppingAction::UserSteppingAction(const G4Step* fStep) 
{
  G4Track* fTrack = fStep->GetTrack();
  G4int StepNo = fTrack->GetCurrentStepNumber();
  if(StepNo >= 10000) fTrack->SetTrackStatus(fStopAndKill);

#ifdef G4ANALYSIS_USE 
#define COLLECT
#endif
#ifdef G4ANALYSIS_USE_ROOT
#ifndef COLLECT
#define COLLECT
#endif
#endif

#ifdef COLLECT
  if (StepNo == 1) {
    if ( (fTrack->GetDefinition()->GetParticleType() == "nucleus") && 
	 !( fTrack->GetDefinition()->GetPDGStable()) && 
	 fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Target" ) {
 //     G4String particleName = fTrack->GetDefinition()->GetParticleName();
      G4double particleName = G4double(fTrack->GetDefinition()->GetPDGEncoding());
      G4double time = fStep->GetPreStepPoint()->GetGlobalTime() ;
	// - fStep->GetPreStepPoint()->GetLocalTime(); // do we have to take out the local time?
      G4double weight =  fStep->GetPreStepPoint()->GetWeight() ;
      // get the analysis manager
      exrdmAnalysisManager* analysis = exrdmAnalysisManager::getInstance();
      //      G4cout << particleName << " " << weight << " " << time/s << G4endl;
      analysis->AddIsotope(particleName, weight, time);
    }
    if (fTrack->GetTrackID() != 1 ){
      if (fTrack->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay") {
	// emitted particles
	if (fTrack->GetDefinition()->GetParticleType() != "nucleus") {
	//  G4String particleName = fTrack->GetDefinition()->GetParticleName();
	  G4double particleName = G4double (fTrack->GetDefinition()->GetPDGEncoding());
	  G4double time = fStep->GetPreStepPoint()->GetGlobalTime() ; 
	  //- fStep->GetPreStepPoint()->GetLocalTime();
	  G4double weight = fStep->GetPreStepPoint()->GetWeight() ;   
	  G4double energy = fStep->GetPreStepPoint()->GetKineticEnergy();
	  //
	  exrdmAnalysisManager::getInstance()->AddParticle(particleName, energy, weight, time);
	}
      }
    }
  }
  // energy deposition
  if (fTrack->GetTrackID() != 1 ) {
    //    if (fTrack->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay") {	
      if (fStep->GetTotalEnergyDeposit() ) {
	G4double time = fStep->GetPreStepPoint()->GetGlobalTime() ;
	// - fStep->GetPreStepPoint()->GetLocalTime(); 
	//    G4double time = fStep->GetPreStepPoint()->GetLocalTime(); // time since the track was created
	G4double edep = fStep->GetTotalEnergyDeposit();
	G4double weight = fStep->GetPreStepPoint()->GetWeight() ; 
	if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Detector") edep = -edep;
	exrdmAnalysisManager::getInstance()->AddEnergy(edep,weight,time);
      }
      //    }
  }
#endif 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


