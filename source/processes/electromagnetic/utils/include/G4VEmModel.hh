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
// GEANT4 Class header file
//
//
// File name:     G4VEmModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.01.2002
//
// Modifications: 

//
// Class Description: 
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#ifndef G4VEmModel_h
#define G4VEmModel_h 1

#include "globals.hh"
#include "g4std/vector"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

class G4VEmModel 
{

public:

  G4VEmModel() {};

  virtual ~G4VEmModel() {};

  virtual G4double HighEnergyLimit(const G4ParticleDefinition*,
                                   const G4Material*) = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition*,
                                  const G4Material*) = 0;

  virtual void SetHighEnergyLimit(const G4Material*, G4double) = 0;
 
  virtual void SetLowEnergyLimit(const G4Material*, G4double) = 0;

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4Material*) = 0;
 
  virtual G4bool IsInCharge(const G4ParticleDefinition*,
			    const G4Material*) = 0;

  virtual G4double ComputeDEDX(const G4Material*,
                               const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     G4double cutEnergy) = 0;

  virtual G4double CrossSection(const G4Material*,
                                const G4ParticleDefinition*,
                                      G4double kineticEnergy,
                                      G4double cutEnergy,
                                      G4double maxEnergy) = 0;

  virtual G4std::vector<G4DynamicParticle*>* SampleSecondary(
                                const G4Material*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double tmax) = 0;

  virtual G4std::vector<G4DynamicParticle*>* DeexciteMedium(
                                const G4Material*,
                                const G4DynamicParticle*,
                                      G4double energyLoss,
                                      G4double cutEnergy) {return 0;};

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle) 
                                {return dynParticle->GetKineticEnergy();};

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                            G4double kineticEnergy) 
                                {return kineticEnergy;};

private:

  // hide assignment operator 
     G4VEmModel & operator=(const  G4VEmModel &right);
     G4VEmModel(const  G4VEmModel&);

};

#endif




