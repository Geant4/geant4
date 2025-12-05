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
#ifndef G4QuasiOpticalData_h
#define G4QuasiOpticalData_h

#include "G4ThreeVector.hh"

// Common data structure used for Cerenkov and Scintillation photon sampling
struct G4QuasiOpticalData
{
  std::size_t mat_index{};  // Index of material in the material table
  G4int num_photons{};  // Number of optical photons to generate
  G4double charge{};  // Charge of the parent particle
  G4double step_length{};  // Step length of the parent track
  G4double pre_velocity{};  // Velocity at the pre-step point
  G4double delta_velocity{};  // Change in velocity between pre- and post-step
  G4ThreeVector delta_position{};  // Displacement vector over the step
};

#endif  // G4QuasiOpticalData_h
