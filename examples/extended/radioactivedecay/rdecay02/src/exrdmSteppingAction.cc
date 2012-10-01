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
/// \file radioactivedecay/rdecay02/src/exrdmSteppingAction.cc
/// \brief Implementation of the exrdmSteppingAction class
//
#include "G4ios.hh"

#include "exrdmSteppingAction.hh"
#include "exrdmAnalysisManager.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"
#include "G4ParticleDefinition.hh"

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
  G4double time = fStep->GetPreStepPoint()->GetGlobalTime();
  G4double weight =  fStep->GetPreStepPoint()->GetWeight();
  G4double edep = fStep->GetTotalEnergyDeposit();
  exrdmAnalysisManager* analysis = exrdmAnalysisManager::GetInstance();
  if (StepNo == 1) { //beginning of step
        G4double pid=G4double(fTrack->GetDefinition()->GetPDGEncoding());
        G4double energy = fStep->GetPreStepPoint()->GetKineticEnergy();
        G4ParticleDefinition* thePartDef = fTrack->GetDefinition();
        G4String partType= fTrack->GetDefinition()->GetParticleType();
        if (( partType == "nucleus") &&  !(thePartDef->GetPDGStable()) &&
       fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Target" ) {
      analysis->AddIsotope(pid, weight, time);
    }
    if (fTrack->GetTrackID() != 1 ){
      //Radioactive decay products
      if (fTrack->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay") {
        //all products
            G4int Z=thePartDef->GetAtomicNumber();
            G4int A=thePartDef->GetAtomicMass();
            analysis->AddDecayProduct(pid, Z,  A, energy, time, weight);
        // emitted particles except nuclei
            if ( partType!= "nucleus") {
                    analysis->AddParticle(pid, energy, weight, time);
            }
      }
    }
  }
  // energy deposition
  if (fTrack->GetTrackID() != 1  && edep) {
   if (fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "Detector")
     edep = -edep;
   analysis->AddEnergy(edep,weight,time);
  }
#endif 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


