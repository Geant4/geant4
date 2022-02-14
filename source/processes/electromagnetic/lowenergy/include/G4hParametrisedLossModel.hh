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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4hParametrisedLossModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 20 July 2000
//
// Modifications: 
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Low energy hadrons/ions ionisation parameterisation
// Further documentation available from http://www.ge.infn.it/geant4/lowE/

// -------------------------------------------------------------------
//

#ifndef G4hParametrisedLossModel_h
#define G4hParametrisedLossModel_h 1

#include "globals.hh"
#include "G4VLowEnergyModel.hh"
#include "G4VhElectronicStoppingPower.hh"

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Material;

class G4hParametrisedLossModel : public G4VLowEnergyModel
{
public:
  explicit G4hParametrisedLossModel(const G4String& name) ;

  ~G4hParametrisedLossModel() ;

  G4double TheValue(const G4DynamicParticle* particle,
	            const G4Material* material) override;

  G4double TheValue(const G4ParticleDefinition* aParticle,
	            const G4Material* material,
                          G4double kineticEnergy) override;

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                           const G4Material* material) const override;

  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                          const G4Material* material) const override;

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle) const override;

  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) const override;

  G4bool IsInCharge(const G4DynamicParticle* particle,
		    const G4Material* material) const override;

  G4bool IsInCharge(const G4ParticleDefinition* aParticle,
		    const G4Material* material) const override;

  G4String ModelName() const {return modelName;};

  G4hParametrisedLossModel(G4hParametrisedLossModel &) = delete;
  G4hParametrisedLossModel & operator=(const G4hParametrisedLossModel &right) = delete;

private:
  void InitializeMe();

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy);

  G4bool MolecIsInZiegler1988(const G4Material* material) ;

  void SetExpStopPower125(G4double value) {expStopPower125 = value;};

  G4double ChemicalFactor(G4double kineticEnergy, G4double eloss125) const;

  // Pointer to the parametrisation class
  G4VhElectronicStoppingPower* eStopingPowerTable;
  G4String modelName;

  G4double theZieglerFactor; // Factor to convert the Stopping Power
                             // unit [ev/(10^15 atoms/cm^2]
                             // into the Geant4 dE/dx unit
  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  G4double expStopPower125;        // Experimental Stopping power at 125keV
};

#endif
