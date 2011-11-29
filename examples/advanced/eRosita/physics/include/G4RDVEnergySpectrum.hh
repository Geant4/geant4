//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4RDVEnergySpectrum
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

#ifndef G4RDVENERGYSPECTRUM_HH
#define G4RDVENERGYSPECTRUM_HH 1

#include "globals.hh"

class G4ParticleDefinition;

class G4RDVEnergySpectrum 
{

public:

  G4RDVEnergySpectrum() {};

  virtual ~G4RDVEnergySpectrum() {};

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
  G4RDVEnergySpectrum(const G4RDVEnergySpectrum&);
  G4RDVEnergySpectrum& operator=(const G4RDVEnergySpectrum &right);

};

#endif

