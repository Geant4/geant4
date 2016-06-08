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
// File name:     G4VEnergySpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 29 September 2001
//
// Modifications: 
//
// -------------------------------------------------------------------

// Class Description: 
//
// Abstract interface for the energy spectrum of secondary particles in
// electromagnetic processes. 
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4VENERGYSPECTRUM_HH
#define G4VENERGYSPECTRUM_HH 1

#include "globals.hh"

class G4ParticleDefinition;

class G4VEnergySpectrum 
{

public:

  G4VEnergySpectrum() {};

  virtual ~G4VEnergySpectrum() {};

  virtual G4double Probability(G4int Z,
			       G4double minKineticEnergy,
			       G4double maxKineticEnergy,
                               G4double kineticEnergy,
                               G4int shell = 0,
			       const G4ParticleDefinition* pd = 0) const = 0;

  virtual G4double AverageEnergy(G4int Z,
				 G4double minKineticEnergy,
				 G4double maxKineticEnergy,
				 G4double kineticEnergy,
				 G4int shell = 0,
				 const G4ParticleDefinition* pd = 0) const = 0;

  virtual G4double SampleEnergy(G4int Z,
				G4double minKineticEnergy,
				G4double maxKineticEnergy,
				G4double kineticEnergy,
				G4int shell = 0,
				const G4ParticleDefinition* pd = 0) const = 0;

  virtual G4double MaxEnergyOfSecondaries(G4double kineticEnergy,
					  G4int Z = 0,
					  const G4ParticleDefinition* pd = 0) const = 0;
  
  virtual G4double Excitation(G4int Z, G4double kineticEnergy) const = 0; 

  virtual void PrintData() const = 0;

protected:

private:

  // Hide copy constructor and assignment operator 
  G4VEnergySpectrum(const G4VEnergySpectrum&);
  G4VEnergySpectrum& operator=(const G4VEnergySpectrum &right);

};

#endif

