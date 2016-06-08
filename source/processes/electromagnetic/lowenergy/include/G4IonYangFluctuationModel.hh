//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     IonYangFluctuationModel
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
// Q.Yang et al., NIM B61(1991)149-155.
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// and Physics Reference Manual

// -------------------------------------------------------------------
//

#ifndef G4IonYangFluctuationModel_h
#define G4IonYangFluctuationModel_h 1

#include "G4VLowEnergyModel.hh"

class G4IonYangFluctuationModel : public G4VLowEnergyModel
{

public: 

  G4IonYangFluctuationModel(const G4String& name);

  ~G4IonYangFluctuationModel();

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

  G4double YangFluctuationModel(const G4Material* material, 
                                      G4double kineticEnergy,
                                      G4double particleMass,
                                      G4double charge) const;

};

#endif
