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
// $Id: G4eIonisationElectronSpectrum.hh,v 1.2 2001-10-10 16:45:56 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eIonisationElectronSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 September 2001
//
// Modifications: 
// 10.10.01 MGP             Revision to improve code quality and 
//                          consistency with design
//
// -------------------------------------------------------------------

// Class Description: 
// Provides various integration over delta-electron spectrum for e- 
// ionisation process

// -------------------------------------------------------------------

#ifndef G4EIONISATIONELECTRONSPECTRUM_HH
#define GG4EIONISATIONELECTRONSPECTRUM_HH 1

#include "G4VEnergySpectrum.hh"
 
class G4eIonisationParameters;
class G4DataVector;

class G4eIonisationElectronSpectrum : public G4VEnergySpectrum
{

public:

  G4eIonisationElectronSpectrum();

  ~G4eIonisationElectronSpectrum();

  G4double Probability(G4int Z, G4double tmin, G4double tmax, 
		       G4double kineticEnergy, G4int shell,
		       const G4ParticleDefinition* pd=0) const;

  G4double AverageEnergy(G4int Z, G4double tmin, G4double tmax,
			 G4double kineticEnergy, G4int shell,
			 const G4ParticleDefinition* pd=0) const;

  G4double SampleEnergy(G4int Z, G4double tmin, G4double tmax,
			G4double kineticEnergy, G4int shell,
			const G4ParticleDefinition* pd=0) const;

  G4double MaxEnergyOfSecondaries(G4double kineticEnergy,
                                  G4int Z = 0,
				  const G4ParticleDefinition* pd=0) const
  { return 0.5*kineticEnergy; }
  
  void PrintData() const;

protected:

private:

  G4double IntSpectrum(size_t n, G4double tmin, G4double tmax, G4double b,
                       const G4DataVector& p) const; 

  G4double AverageValue(size_t n, G4double tmin, G4double tmax, G4double b,
                       const G4DataVector& p) const; 

  G4double Function(size_t n, G4double e, G4double b, 
                       const G4DataVector& p) const; 

  G4double MaxFunction(size_t n, G4double tmin, G4double tmax, G4double b, 
                       const G4DataVector& p) const; 

  // Hide copy constructor and assignment operator 
  G4eIonisationElectronSpectrum(const  G4eIonisationElectronSpectrum&);
  G4eIonisationElectronSpectrum & operator = (const G4eIonisationElectronSpectrum &right);

private:

  G4eIonisationParameters* theParam;
  G4double                 lowestE;
  G4int                    verbose;      

};


#endif
