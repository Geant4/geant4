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
// Ion effective charge model
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4hIonEffChargeSquare_h
#define G4hIonEffChargeSquare_h 1

#include "G4VLowEnergyModel.hh"

class G4hIonEffChargeSquare : public G4VLowEnergyModel
{

public:

  G4hIonEffChargeSquare(const G4String& name) ;

  ~G4hIonEffChargeSquare() ;

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

  G4double IonEffChargeSquare(const G4Material* material, 
                                    G4double kineticEnergy,
                                    G4double particleMass,
                                    G4double ionCharge) const;
  // This method returns ion effective charge square parametrised according to
  // J.F.Ziegler, J.P. Biersack, U. Littmark
  // The Stopping and Range of Ions in Matter,
  // Vol.1, Pergamon Press, 1985

  const G4double protonMass;
  const G4double theHeMassAMU;

};

#endif
