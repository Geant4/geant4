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

#ifndef G4LowEnergyDeltaGen_h
#define G4LowEnergyDeltaGen_h 1

#include "G4VEMSecondaryGenerator.hh" 
#include "globals.hh"
#include "G4DataVector.hh"
#include "g4std/map"
#include "g4std/vector"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4ParticleDefinition;
class G4DynamicParticle;
class G4Material;
class G4ParticleChange;
class G4VDataSetAlgorithm;
class G4VEMDataSet;
class G4ShellEMDataSet;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4LowEnergyDeltaGen : public G4VEMSecondaryGenerator
{

public:

  G4LowEnergyDeltaGen(const G4String& name = "LowEnDeltaGen");

  virtual ~G4LowEnergyDeltaGen();

  void Initialize();

  void Clear();

  G4double Probability(G4int, G4int, G4double, G4double, G4double) const;

  G4double AverageEnergy(G4int, G4int, G4double, G4double) const;

  G4double CrossSectionWithCut(G4int,G4double,G4double,G4double) const 
           {return 0.0;};

  void GenerateSecondary(const G4DynamicParticle*, G4ParticleChange*,
                               G4int, G4int, G4double, G4double);

  G4double MinSecondaryEnergy(const G4Material*) const 
           {return 0.0;};

  //  void SetMinSecondaryEnergy(G4double val) {lowestEnergyGamma = val;};

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double kineticEnergy) const
           {return 0.5*kineticEnergy;};

  void SetHighEnergyLimit(G4double val) {highEnergyLimit = val;};

  G4double HighEnergyLimit(const G4ParticleDefinition* aParticle,
                           const G4Material*) const
           {return highEnergyLimit;};

  void SetLowEnergyLimit(G4double val) {lowEnergyLimit = val;};
 
  G4double LowEnergyLimit(const G4ParticleDefinition*,
                          const G4Material*) const
           {return lowEnergyLimit;};

  G4double HighEnergyLimit(const G4ParticleDefinition*) const
           {return highEnergyLimit;};

  G4double LowEnergyLimit(const G4ParticleDefinition*) const
           {return lowEnergyLimit;};
 
  G4bool IsInCharge(const G4ParticleDefinition* aParticle,
		    const G4Material*) const
           {return (aParticle == G4Electron::Electron());};

  G4bool IsInCharge(const G4ParticleDefinition* aParticle) const
           {return (aParticle == G4Electron::Electron());};

  virtual void PrintData() const;

protected:

private:

  // hide assignment operator 
  G4LowEnergyDeltaGen(const  G4LowEnergyDeltaGen&);
  G4LowEnergyDeltaGen & operator =
                              (const  G4LowEnergyDeltaGen &right);

  const G4DataVector* FindParameters(G4int, G4int, G4double) const;

  G4double FindEnergyLimit1(G4int, G4int) const;

  G4double FindEnergyLimit2(G4int, G4int) const;

  G4double IntSpectrum1(size_t, size_t, G4double, G4double, G4double,
                        const G4DataVector&) const;

  G4double IntSpectrum2(size_t, G4double, G4double, const G4DataVector&) const;

  G4double IntSpectrum3(size_t, size_t, G4double, G4double, G4double) const;

  G4double Spectrum1(size_t, G4double, G4double, const G4DataVector&) const;

  G4double Spectrum2(G4double, const G4DataVector&) const;

  // Parameters of the energy spectra
  G4DataVector activeZ;

  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyLimit1;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyLimit2;
  G4std::map<G4int,G4VEMDataSet*,G4std::less<G4int> > param;

  size_t length_a;
  size_t length_p;
  size_t length_z;

  // The interpolation algorithm
  const G4VDataSetAlgorithm* interpolation;
  const G4VDataSetAlgorithm* interpolation1;

  const G4ShellEMDataSet* bindingData;

  // Limit of the model validity
  G4double lowEnergyLimit;    
  G4double highEnergyLimit;    

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif



