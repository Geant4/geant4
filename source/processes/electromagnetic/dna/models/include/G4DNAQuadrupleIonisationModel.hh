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
// G4DNAQuadrupleIonisationModel.hh
//
//  Created at 2024/04/03 (Thu.)
//  Author: Shogo OKADA @KEK-CRC (shogo.okada@kek.jp)
//
//  Reference: J.Meesungnoen et. al, DOI: 10.1021/jp058037z
//

#ifndef G4DNA_QUADRUPLE_IONISATION_MODEL_HH_
#define G4DNA_QUADRUPLE_IONISATION_MODEL_HH_

#include "G4DNADoubleIonisationModel.hh"
#include "G4VEmModel.hh"

class G4DNAQuadrupleIonisationModel : public G4DNADoubleIonisationModel {
public:
  // constructor
  G4DNAQuadrupleIonisationModel(
    const G4ParticleDefinition* p = nullptr,
    const G4String& model_name = "G4DNAQuadrupleIonisationModel");

  // destructor
  ~G4DNAQuadrupleIonisationModel() override = default;

  void Initialise(const G4ParticleDefinition* particle,
                  const G4DataVector&) override;

  G4double CrossSectionPerVolume(
    const G4Material* material, const G4ParticleDefinition* pdef,
    G4double ekin, G4double, G4double) override;

  void SampleSecondaries(
    std::vector<G4DynamicParticle*>* vsec, const G4MaterialCutsCouple* couple,
    const G4DynamicParticle* particle, G4double, G4double) override;
};

#endif // G4DNA_QUADRUPLE_IONISATION_MODEL_HH_
