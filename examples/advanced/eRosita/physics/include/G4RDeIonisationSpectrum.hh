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
// File name:     G4RDeIonisationSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 September 2001
//
// Modifications: 
// 10.10.01 MGP             Revision to improve code quality and 
//                          consistency with design
// 29.11.01  V.Ivanchenko   Parametrisation is updated
// 30.05.02  VI             Add inline functions
//
// -------------------------------------------------------------------

// Class Description: 
// Provides various integration over delta-electron spectrum for e- 
// ionisation process. Spectrum is parametrised accourding to 
// EEDL database, details are described in the Physics Reference Manual
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDEIONISATIONSPECTRUM_HH
#define G4RDEIONISATIONSPECTRUM_HH 1

#include "G4RDVEnergySpectrum.hh"
#include "G4RDeIonisationParameters.hh"
 
//class G4RDeIonisationParameters;
class G4DataVector;

class G4RDeIonisationSpectrum : public G4RDVEnergySpectrum
{

public:

  G4RDeIonisationSpectrum();

  ~G4RDeIonisationSpectrum();

  G4double Probability(G4int Z, G4double tMin, G4double tMax, 
		       G4double kineticEnergy, G4int shell,
		       const G4ParticleDefinition* pd=0) const;

  G4double AverageEnergy(G4int Z, G4double tMin, G4double tMax,
			 G4double kineticEnergy, G4int shell,
			 const G4ParticleDefinition* pd=0) const;

  G4double SampleEnergy(G4int Z, G4double tMin, G4double tMax,
			G4double kineticEnergy, G4int shell,
			const G4ParticleDefinition* pd=0) const;

  G4double MaxEnergyOfSecondaries(G4double kineticEnergy,
                                  G4int Z = 0,
				  const G4ParticleDefinition* pd=0) const;

  G4double Excitation(G4int Z, G4double e) const;
  
  void PrintData() const;

protected:

private:

  G4double IntSpectrum(G4double xMin, G4double xMax,
                         const G4DataVector& p) const; 
  
  G4double AverageValue(G4double xMin, G4double xMax,
			const G4DataVector& p) const; 
  
  G4double Function(G4double x, const G4DataVector& p) const; 
  
  // Hide copy constructor and assignment operator 
  G4RDeIonisationSpectrum(const  G4RDeIonisationSpectrum&);
  G4RDeIonisationSpectrum & operator = (const G4RDeIonisationSpectrum &right);
  
private:
  
  G4RDeIonisationParameters* theParam;
  G4double                 lowestE;
  G4double                 factor;
  G4int                    iMax;            
  G4int                    verbose;
};


inline G4double G4RDeIonisationSpectrum::Function(G4double x, 
				          const G4DataVector& p) const
{
  G4double f = 1.0 - p[0] - p[iMax]*x 
             + x*x*(1.0 - p[iMax] + (1.0/(1.0 - x) - p[iMax])/(1.0 - x))
             + 0.5*p[0]/x;

  return f;
} 

inline G4double G4RDeIonisationSpectrum::Excitation(G4int Z, G4double e) const 
{
  return theParam->Excitation(Z, e);
}


#endif
