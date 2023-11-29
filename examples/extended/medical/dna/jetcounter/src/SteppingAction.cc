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
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction *event)
  : G4UserSteppingAction()
  , fEventAction(event) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step *step) {
  // /////////////////////////////////////////
  // killing primary particles (any heavy charged particle) that hit the
  // collimator
  auto part_name = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "collPhys" &&
      step->GetTrack()->GetParticleDefinition()->GetParticleName() == "alpha") {
    step->GetTrack()->SetTrackStatus(fStopAndKill);
    fEventAction->InvalidEvent();
    return; // nothing more to do in such case
  }

  // /////////////////////////////////////////
  // recording energy of the primary particle at three different locations:
  // 1. Boundary between World and Mylar wall, i.e. initial energy of the
  // particle.
  // 2. Boundary between World and silicon detector, i.e. final energy of the
  // particle.
  // 3. Boundary between wall inner layer (innerWall) and Target, i.e.
  // interaction energy of the particle.

  if (step->GetTrack()->GetParticleDefinition()->GetParticleName() == "alpha") {
    // 1. and 2. - particle crossing the boundary between World volume
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() ==
        "worldPhys") {
      auto physVol = step->GetPostStepPoint()->GetPhysicalVolume();
      // the above may be a null pointer as the particle can leave world, so it
      // has to be tested:
      if (physVol) {
        // 1. to Mylar wall
        if (physVol->GetName() == "wallPhys") {
          fEventAction->RecordInitialEnergy(
              step->GetTrack()->GetKineticEnergy());
        }
        // 2. to silicon detector
        if (physVol->GetName() == "enDetPhys") {
          fEventAction->RecordFinalEnergy(step->GetTrack()->GetKineticEnergy());
        }
      }
    }

    // 3. - particle crossing the boundary from wall inner layer to the target
    if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() ==
        "innerWallPhys") {
      if (step->GetPostStepPoint()->GetPhysicalVolume()->GetName() ==
          "targetPhys") {
        fEventAction->RecordInteractionEnergy(
            step->GetTrack()->GetKineticEnergy());
      }
    }
  }

  // ////////////////////////////////////////
  // counting individual ionization acts
  // only ionisations in the target will be recorded
  if (step->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "targetPhys") {
    if (step->GetTotalEnergyDeposit() != 0) {
      G4String process_name =
          step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

      // counting ionisations caused only by electrons and alpha particles:
      //            if(process_name=="e-_G4DNAIonisation" ||
      //            process_name=="alpha_G4DNAIonisation");
      //                fEventAction->AddIonisation();

      // more general approach - counting ionisations caused by any particle
      // (including also protons, GenericIons, etc.):
      auto pos = process_name.find(
          "Ionisation"); // looking for "Ionisation" substring in process_name
      if (pos != std::string::npos)    // found!
        fEventAction->AddIonisation(); // increase cluster size
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
