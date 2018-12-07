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
// File name:     G4EBremsstrahlungSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
//
// Creation date: 27 September 2001
//
// Modifications:
// 10.10.01  MGP  Revision to improve code quality and consistency with design
// 29.11.01  V.Ivanchenko    Parametrisation is updated
// 21.02.03  V.Ivanchenko    Energy bins are defined in the constructor
// 28.02.03  V.Ivanchenko    Filename is defined in the constructor
// 25.05.03  MGP             Data member xp contained by value, not a reference
//
// -------------------------------------------------------------------

// Class Description:
// Provides various integration over gamma spectrum of e- Bremsstrahlung.
// Parametrisation is described in Physics Reference Manual based on
// data from EEDL database.
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4RDEBREMSSTRAHLUNGSPECTRUM_HH
#define G4RDEBREMSSTRAHLUNGSPECTRUM_HH 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4RDVEnergySpectrum.hh"

class G4RDBremsstrahlungParameters;

class G4RDeBremsstrahlungSpectrum : public G4RDVEnergySpectrum
{
public:

  G4RDeBremsstrahlungSpectrum(const G4DataVector& bins,const G4String& name);

  ~G4RDeBremsstrahlungSpectrum();

  G4double Probability(G4int Z,
                       G4double tMin,
                       G4double tMax,
                       G4double kineticEnergy,
                       G4int shell=0,
		       const G4ParticleDefinition* pd=0) const;

  G4double AverageEnergy(G4int Z,
                         G4double tMin,
                         G4double tMax,
                         G4double kineticEnergy,
                         G4int shell=0,
			 const G4ParticleDefinition* pd=0) const;

  G4double SampleEnergy(G4int Z,
                        G4double tMin,
                        G4double tMax,
                        G4double kineticEnergy,
                        G4int shell=0,
			const G4ParticleDefinition* pd=0) const;

  G4double MaxEnergyOfSecondaries(G4double kineticEnergy,
                                  G4int Z = 0,
				  const G4ParticleDefinition* pd=0) const;

  G4double Excitation(G4int Z, G4double kineticEnergy) const;

  void PrintData() const;

private:

  G4double IntSpectrum(G4double xMin, G4double xMax,
                         const G4DataVector& p) const;

  G4double AverageValue(G4double xMin, G4double xMax,
			const G4DataVector& p) const;

  G4double Function(G4double x, const G4DataVector& p) const;


  // Hide copy constructor and assignment operator
  G4RDeBremsstrahlungSpectrum(const  G4RDeBremsstrahlungSpectrum&);
  G4RDeBremsstrahlungSpectrum & operator = (const G4RDeBremsstrahlungSpectrum &right);

  G4RDBremsstrahlungParameters* theBRparam;
  G4double                    lowestE;
  size_t                      length;
  G4int                       verbose;
  // const G4DataVector&         xp;
  const G4DataVector          xp;
};

#endif
