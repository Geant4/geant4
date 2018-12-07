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
//
/// \file field/field04/src/F04StackingAction.cc
/// \brief Implementation of the F04StackingAction class
//

#include "F04StackingAction.hh"

#include "G4RunManager.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04StackingAction::F04StackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04StackingAction::~F04StackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
      F04StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ParticleDefinition* particleType = aTrack->GetDefinition();
  const G4String name = particleType->GetParticleName();

  // keep primary particle
  if (aTrack->GetParentID() == 0) return fUrgent;

  if (particleType != G4Proton::ProtonDefinition()     &&
      particleType != G4Neutron::NeutronDefinition()   &&
      particleType != G4KaonPlus::KaonPlusDefinition() &&
      particleType != G4PionPlus::PionPlusDefinition() &&
      particleType != G4MuonPlus::MuonPlusDefinition() &&
      particleType != G4Positron::PositronDefinition()) return fKill;

  if (name != "pi+" && name != "mu+") return fKill;

//  if (name == "mu+") G4RunManager::GetRunManager()->rndmSaveThisEvent();

  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04StackingAction::NewStage() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04StackingAction::PrepareNewEvent() {}
