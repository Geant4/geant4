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
// $Id: G4eBremstrahlungElectronSpectrum.hh,v 1.1 2001-10-10 11:51:01 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// 10 Oct 2001  MGP  Revision to improve code quality and consistency with design
//
// -------------------------------------------------------------------

// Class Description: 
// Provides various integration over gamma spectrum of e- Bremsstrahlung  

// -------------------------------------------------------------------

#ifndef G4EBREMSSTRAHLUNGSPECTRUM_HH
#define G4EBREMSSTRAHLUNGSPECTRUM_HH 1

#include "globals.hh"
#include "G4VEnergySpectrum.hh" 

class G4BremsstrahlungParameters;

class G4eBremsstrahlungSpectrum : public G4VEnergySpectrum
{
public:

  G4eBremsstrahlungSpectrum();

  ~G4eBremsstrahlungSpectrum();

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
				  const G4ParticleDefinition* pd=0) const
  { return kineticEnergy; }

  void PrintData() const;

private:

  // Hide copy constructor and assignment operator 
  G4eBremsstrahlungSpectrum(const  G4eBremsstrahlungSpectrum&);
  G4eBremsstrahlungSpectrum & operator = (const G4eBremsstrahlungSpectrum &right);

  G4BremsstrahlungParameters* theBRparam;
  G4double                    lowestE;

};

#endif
