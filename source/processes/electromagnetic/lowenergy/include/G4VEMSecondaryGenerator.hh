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
// File name:     G4VEMSecondaryGenerator
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 August 2001
//
// Modifications: 
//
// Class Description: 
//
// Abstract class provided the interface to generators of secondary
// particles for electromagnetic processes. 

// -------------------------------------------------------------------
//

#ifndef G4VEMSecondaryGenerator_h
#define G4VEMSecondaryGenerator_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;
class G4ParticleChange;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4VEMSecondaryGenerator 
{

public:

  G4VEMSecondaryGenerator(const G4String& name) {generatorName = name;};

  virtual ~G4VEMSecondaryGenerator() {};

  virtual void Initialize() = 0;

  virtual void Clear() = 0;

  virtual G4double Probability(G4int atomicNumber,
			       G4double kineticEnergy,
                               G4double tmin,
                               G4double tmax) const = 0;

  virtual G4double AverageEnergy(G4int atomicNumber,
			       G4double kineticEnergy,
                               G4double tcut) const = 0;

  virtual G4double CrossSectionWithCut(G4int atomicNumber,
			       G4double kineticEnergy,
                               G4double tmin,
                               G4double tmax) const = 0;

  virtual void GenerateSecondary(const G4DynamicParticle* aParticle,
	                               G4ParticleChange* theChange,
                                       G4int atomicNumber,
                                       G4int shellNumber,
                                       G4double tmin,
                                       G4double tmax) = 0;

  virtual G4double MinSecondaryEnergy(const G4Material* material) const = 0;

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition* aParticle,
		                      const G4Material* material,
                                            G4double kineticEnergy) const = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                                   const G4Material* material) const = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition* aParticle,
                                  const G4Material* material) const = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition* aParticle)
                                   const = 0;
 
  virtual G4double LowEnergyLimit(const G4ParticleDefinition* aParticle) 
                                  const = 0;
 
  virtual G4bool IsInCharge(const G4ParticleDefinition* aParticle,
			    const G4Material* material) const = 0;

  virtual G4bool IsInCharge(const G4ParticleDefinition* aParticle)
			    const = 0;

  virtual void PrintData() const = 0;

  G4String GeneratorName() const {return generatorName;};

  void SetVerbose(G4int val) {verbose = val;};

protected:

  G4String generatorName;
  G4int    verbose;

private:

  // hide assignment operator 
  G4VEMSecondaryGenerator(const  G4VEMSecondaryGenerator&);
  G4VEMSecondaryGenerator & operator=(const  G4VEMSecondaryGenerator &right);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

