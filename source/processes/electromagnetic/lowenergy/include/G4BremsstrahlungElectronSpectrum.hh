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
// File name:     G4BremsstrahlungElectronSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 September 2001
//
// Modifications: 
//
// -------------------------------------------------------------------

// Class Description: 
//
// Provides various integration over gamma spectrum of e- bremsstrahlung  
//
// Class Description: End 

// -------------------------------------------------------------------
//

#ifndef G4BremsstrahlungElectronSpectrum_h
#define G4BremsstrahlungElectronSpectrum_h 1

#include "G4VEnergySpectrum.hh" 
#include "G4BremsstrahlungParameters.hh" 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4BremsstrahlungElectronSpectrum : public G4VEnergySpectrum
{

public:

  G4BremsstrahlungElectronSpectrum();

  ~G4BremsstrahlungElectronSpectrum();

  G4double Probability(G4int Z, 
                       G4double tmin, 
                       G4double tmax, 
                       G4double kineticEnergy, 
                       G4int shell=0, 
                 const G4ParticleDefinition* pd=0) const;

  G4double AverageEnergy(G4int Z, 
                         G4double tmin, 
                         G4double tmax,
                         G4double kineticEnergy,
                         G4int shell=0, 
                   const G4ParticleDefinition* pd=0) const;

  G4double SampleEnergy(G4int Z, 
                        G4double tmin, 
                        G4double tmax,
                        G4double kineticEnergy,
                        G4int shell=0, 
                  const G4ParticleDefinition* pd=0) const;

  G4double MaxEnergyOfSecondaries(G4double kineticEnergy,
                                  G4int Z = 0,
                            const G4ParticleDefinition* pd=0) const
                       {return kineticEnergy;};

  void PrintData() const {theBRparam->PrintData();};

protected:

private:

  // hide assignment operator 
  G4BremsstrahlungElectronSpectrum(const  G4BremsstrahlungElectronSpectrum&);
  G4BremsstrahlungElectronSpectrum & operator =
                             (const G4BremsstrahlungElectronSpectrum &right);

private:

  G4BremsstrahlungParameters* theBRparam;
  G4double                    lowestE;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif



