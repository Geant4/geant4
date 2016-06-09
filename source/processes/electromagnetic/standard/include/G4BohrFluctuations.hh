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
// $Id: G4BohrFluctuations.hh,v 1.1 2004/12/01 17:36:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BohrFluctuations
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 02.04.2003
//
// Modifications: 
//
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
//
// Class Description:
//
// Implementation of Gaussion energy loss Fluctuations

// -------------------------------------------------------------------
//

#ifndef G4BohrFluctuations_h
#define G4BohrFluctuations_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"

class G4BohrFluctuations : public G4VEmFluctuationModel
{

public:

  G4BohrFluctuations(const G4String& nam = "BohrFluc");

  virtual ~G4BohrFluctuations();

  G4double SampleFluctuations(const G4Material*,
                              const G4DynamicParticle*,
 				    G4double&,
                                    G4double&,
                                    G4double&);

  G4double Dispersion(    const G4Material*,
                          const G4DynamicParticle*,
 				G4double&,
                                G4double&);

  void InitialiseMe(const G4ParticleDefinition*);

protected:

private:

  // hide assignment operator
  G4BohrFluctuations & operator=(const  G4BohrFluctuations &right);
  G4BohrFluctuations(const  G4BohrFluctuations&);

  const G4ParticleDefinition* particle;

  G4double particleMass;
  G4double chargeSquare;

  G4double minNumberInteractionsBohr;
  G4double minFraction;
  G4double xmin;
  G4double minLoss;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline G4double G4BohrFluctuations::Dispersion(
                          const G4Material* material,
                          const G4DynamicParticle* dp,
 				G4double& tmax,
			        G4double& length)
{
  if(!particle) InitialiseMe(dp->GetDefinition());

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

