// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// For information related to this code contact:
// Geant4 Collaboration
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
// Low energy hadrons/ions ionisation parametrisation
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4hParametrisedLossModel_h
#define G4hParametrisedLossModel_h 1

#include "G4VLowEnergyModel.hh"
#include "G4VhElectronicStoppingPower.hh"

class G4hParametrisedLossModel : public G4VLowEnergyModel
{

public:

  G4hParametrisedLossModel(const G4String& name) ;

  ~G4hParametrisedLossModel() ;

  G4double TheValue(const G4DynamicParticle* particle,
	            const G4Material* material);

  G4double TheValue(const G4ParticleDefinition* aParticle,
	            const G4Material* material,
                          G4double kineticEnergy);

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                           const G4Material* material) const;
 
  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                          const G4Material* material) const;

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle) const;
 
  G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) const;
 
  G4bool IsInCharge(const G4DynamicParticle* particle,
		    const G4Material* material) const;

  G4bool IsInCharge(const G4ParticleDefinition* aParticle,
		    const G4Material* material) const;

protected:

private:

  // hide  assignment operator 
  G4hParametrisedLossModel(G4hParametrisedLossModel &);
  G4hParametrisedLossModel & operator=(const G4hParametrisedLossModel &right);

  void InitializeMe();

  G4double StoppingPower(const G4Material* material,
                               G4double kineticEnergy);

  G4bool MolecIsInZiegler1988(const G4Material* material) ;

  void SetExpStopPower125(G4double value) {expStopPower125 = value;};

  G4double ChemicalFactor(G4double kineticEnergy, G4double eloss125) const;

  // Pointer to the parametrisation class
  G4VhElectronicStoppingPower* eStopingPowerTable;

  G4double theZieglerFactor; // Factor to convert the Stopping Power 
                             // unit [ev/(10^15 atoms/cm^2]
                             // into the Geant4 dE/dx unit
  G4String modelName;

  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  
  G4double expStopPower125;        // Experimental Stopping power at 125keV

};

#endif
