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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <fstream>

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction *detector)
    : G4VUserPrimaryGeneratorAction(), fDetector(detector){
  fpParticleGun = std::make_unique<G4ParticleGun>();
  fpParticleGun->SetParticleEnergy(
      -1); // default value - can be overridden in the macro file
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  // user can set the energy in the macro file
  G4double energy = fpParticleGun->GetParticleEnergy();

  // only first check is important as only user can set energy to -1,
  // hopefully...
  if (energy == -1) { // if energy is larger than zero, then the beam is
                      // mono-energetic if the energy is set to zero, then the
                      // energy is randomized, based on spectrum file
    fMonoEnergetic = false;
  }

  if (!fMonoEnergetic) {
    energy = GenerateParticleEnergy();
  }

  fpParticleGun->SetParticleEnergy(energy);
  fpParticleGun->SetParticlePosition(
      GenerateParticlePosition()); // point of emission
  fpParticleGun->SetParticleMomentumDirection(
      GenerateParticleDirection()); // direction of emission

  fpParticleGun->GeneratePrimaryVertex(anEvent); // sending the particle
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::GenerateParticleEnergy() {
  if (fEnergySpectrum_length ==
      0) { // reading the spectrum file if not loaded yet
    std::ifstream fin(fEnergySpectrumFilename);
    G4String len, gain, offset, counts;
    fin >> len >> gain >> offset;
    fEnergySpectrum_length = std::stoi(len);
    fEnergySpectrum_gain = std::stod(gain);
    fEnergySpectrum_offset = std::stod(offset);

    fEnergySpectrum_counts.resize(fEnergySpectrum_length);
    for (G4int i = 0; i < fEnergySpectrum_length; i++) {
      fin >> counts;
      fEnergySpectrum_counts[i] = std::stoi(counts);
    }
  }

  auto max_en_counter = fEnergySpectrum_counts[fEnergySpectrum_length - 1] - 1;
  auto rand_count = 1 + G4int(G4UniformRand() * max_en_counter);

  // primitive search for lowest en_id that gives higher count value than
  // rand_count:
  G4int en_id = 0;
  for (; rand_count > fEnergySpectrum_counts[en_id]; en_id++)
    ;

  G4int left = fEnergySpectrum_counts[en_id - 1];
  G4int right = fEnergySpectrum_counts[en_id];

  G4double en_left =
      (en_id - 1) * fEnergySpectrum_gain + fEnergySpectrum_offset;
  G4double en_right = en_id * fEnergySpectrum_gain + fEnergySpectrum_offset;

  // linear interpolation:
  G4double slope = (en_right - en_left) / (right - left);
  return (en_left + slope * (rand_count - left)) * MeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::GenerateParticlePosition() {
  // the source is an infinitely thin disk of radius r
  G4double r = std::sqrt(G4UniformRand()) * fDetector->GetCollDiameter() / 2.;
  G4double phi = G4UniformRand() * twopi;

  G4double x = fDetector->GetCollExitPosistion() - fDetector->GetCollLength();
  G4double y = r * std::cos(phi);
  G4double z = r * std::sin(phi);

  return {x, y, z};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::GenerateParticleDirection() {
  G4double phi, theta;
  G4double px, py, pz;

  phi = G4UniformRand() * twopi;

  // For convenience theta is measured along the beam axis, which is x-axis.
  // To reduce the number of events, when the projectile hits the collimator,
  // only forward angles in the range [0, max_theta] are considered.
  // If theta is larger, then the projectile will hit the collimator anyway.
  // Even with this restriction only 1 in 4 projectiles pass through the
  // collimator.
  G4double cos_max_theta =
      std::cos(std::atan(fDetector->GetCollDiameter() / fDetector->GetCollLength()));
  theta = std::acos(cos_max_theta + G4UniformRand() * (1 - cos_max_theta));

  px = std::cos(theta);
  py = std::sin(theta) * std::sin(phi);
  pz = std::sin(theta) * std::cos(phi);

  return {px, py, pz};
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
