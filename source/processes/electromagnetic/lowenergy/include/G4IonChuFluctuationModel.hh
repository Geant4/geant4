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
// File name:     IonChuFluctuationModel
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 18 August 2000
//
// Modifications: 
// 18/08/2000  V.Ivanchenko First implementation
//
// Class Description: 
//
// The aproximation of additional ion energy loss fluctuations 
// W.K.Chu, In: Ion Beam Handbook for Material Analysis.
// eds. J.W. Mayer and E. Rimini (Academic Press, New York, 1977).
// Q.Yang et al., NIM B61(1991)149-155.
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

#ifndef G4IonChuFluctuationModel_h
#define G4IonChuFluctuationModel_h 1

#include "G4VLowEnergyModel.hh"

class G4IonChuFluctuationModel : public G4VLowEnergyModel
{

public:

  G4IonChuFluctuationModel(const G4String& name) ;

  ~G4IonChuFluctuationModel() ;

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

  G4double ChuFluctuationModel(const G4Material* material, 
                                     G4double kineticEnergy,
                                     G4double particleMass) const;

  const G4double protonMassAMU;

};

#endif
