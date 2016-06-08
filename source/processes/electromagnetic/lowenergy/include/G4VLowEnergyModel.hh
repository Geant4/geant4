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
// File name:     G4VLowEnergyModel
//
// Author:        Maria Grazia Pia (MariaGrazia.Pia@ge.infn.it)
// 
// Creation date: 7 May 2000
//
// Modifications: 
// 22/05/2000  MGP          Version compliant with design
// 20/07/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// Abstract class for Low Energy Electromagnetic models
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4VLowEnergyModel_h
#define G4VLowEnergyModel_h 1

#include "G4ios.hh"
#include "globals.hh"

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;

class G4VLowEnergyModel 
{

public:

  G4VLowEnergyModel(const G4String& name);

  virtual ~G4VLowEnergyModel();

  virtual G4double TheValue(const G4DynamicParticle* particle,
			    const G4Material* material)  = 0;

  virtual G4double TheValue(const G4ParticleDefinition* aParticle,
			    const G4Material* material,
                                  G4double kineticEnergy) = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                                   const G4Material* material) const = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                                  const G4Material* material) const = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition* aParticle)
                                  const = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) 
                                  const = 0;
 
  virtual G4bool IsInCharge(const G4DynamicParticle* particle,
			    const G4Material* material) const = 0;
 
  virtual G4bool IsInCharge(const G4ParticleDefinition* aParticle,
			    const G4Material* material) const = 0;

protected:

private:

  // hide assignment operator 
     G4VLowEnergyModel & operator=(const  G4VLowEnergyModel &right);
     G4VLowEnergyModel(const  G4VLowEnergyModel&);

};

#endif

