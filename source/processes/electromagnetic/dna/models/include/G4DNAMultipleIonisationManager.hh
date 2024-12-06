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
// G4DNAMultipleIonisationManager.hh
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//

#ifndef G4DNA_MULTIPLE_IONISATION_MANAGER_HH_
#define G4DNA_MULTIPLE_IONISATION_MANAGER_HH_

#include "globals.hh"
#include <vector>

class G4Track;

enum MultipleIonisedModification {
  eDoubleIonisedMolecule    = 3,
  eTripleIonisedMolecule    = 4,
  eQuadrupleIonisedMolecule = 5
};

class G4DNAMultipleIonisationManager {
public:
  G4DNAMultipleIonisationManager() = default;

  ~G4DNAMultipleIonisationManager() = default;

  void CreateMultipleIonisedWaterMolecule(
    MultipleIonisedModification mod, G4int* shell_level,
    const G4Track* incoming_track);

  G4bool CheckShellEnergy(
    MultipleIonisedModification mod, G4double* shell_energy);

  void LoadAlphaParam(const G4String& filepath, G4double Z, G4double A);

  G4double GetAlphaParam(G4double energy);

private:
  G4int num_node_{0};
  std::vector<G4double> Eion_;
  std::vector<G4double> alpha_;

};

#endif // G4DNA_MULTIPLE_IONISATION_MANAGER_HH_
