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
// File name:     G4BohrFluctuation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 02.04.2003
//
// Modifications:
//
//
// Class Description:
//
// Implementation of Gaussion energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4BohrFluctuation_h
#define G4BohrFluctuation_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"

class G4BohrFluctuation : public G4VEmFluctuationModel
{

public:

  G4BohrFluctuation(const G4String& nam = "BohrFluc");

  ~G4BohrFluctuation();

  G4double SampleFluctuations(const G4Material*,
                              const G4DynamicParticle*,
 				    G4double&,
                                    G4double&,
                                    G4double&);

  G4double Dispersion(    const G4Material*,
                          const G4DynamicParticle*,
 				G4double&,
                                G4double&);

  void Initialise(const G4ParticleDefinition*);

protected:

private:

  // hide assignment operator
  G4BohrFluctuation & operator=(const  G4BohrFluctuation &right);
  G4BohrFluctuation(const  G4BohrFluctuation&);

  const G4ParticleDefinition* particle;

  G4double particleMass;
  G4double chargeSquare;

  G4double minNumberInteractionsBohr;
  G4double minFraction;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline G4double G4BohrFluctuation::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
 				G4double& tmax,
			        G4double& length)
{
  G4double electronDensity = material->GetElectronDensity();
  kineticEnergy  = dp->GetKineticEnergy();
  G4double gam   = kineticEnergy/particleMass + 1.0;
  beta2 = 1.0 - 1.0/(gam*gam);
  G4double siga  = (1.0/beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length
                 * electronDensity * chargeSquare;

  return siga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

