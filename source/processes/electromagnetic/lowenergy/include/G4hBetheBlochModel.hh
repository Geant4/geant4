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
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4hBetheBlochModel_h
#define G4hBetheBlochModel_h 1

#include "globals.hh"
#include "G4VLowEnergyModel.hh"

class G4Material;
class G4ParticleDefinition;
class G4DynamicParticle;

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
