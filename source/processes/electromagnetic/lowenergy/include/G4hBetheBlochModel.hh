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
// File name:     G4hBetheBlochModel
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
// Bethe-Bloch ionisation model
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4hBetheBlochModel_h
#define G4hBetheBlochModel_h 1

#include "G4VLowEnergyModel.hh"

class G4hBetheBlochModel : public G4VLowEnergyModel
{

public:

  G4hBetheBlochModel(const G4String& name) ;

  ~G4hBetheBlochModel() ;

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

  G4double BetheBlochFormula(const G4Material* material,
                                   G4double kineticEnergy,
                                   G4double particleMass) const;

  // Low energy limit of the model
  G4double lowEnergyLimit;
  G4double highEnergyLimit;

  // constants needed for the energy loss calculation
  
  const G4double twoln10;
  const G4double bg2lim;
  const G4double taulim;    // energy to start to switch off shell corrections

};

#endif
