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
// $Id: G4ionEffectiveCharge.hh,v 1.6 2005/02/26 22:01:20 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ionEffectiveCharge
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.07.2004
//
// Modifications:
//
//
// Class Description:
//
// This class manages the simulation of effective charge of ions 
// in the assumption of equilibrium between ion shelss and media.
// J.F.Ziegler and J.M.Manoyan, The stopping of ions in compaunds,
// Nucl. Inst. & Meth. in Phys. Res. B35 (1988) 215-228.
//

// -------------------------------------------------------------------
//

#ifndef G4ionEffectiveCharge_h
#define G4ionEffectiveCharge_h 1

#include "globals.hh"

class G4Material;
class G4ParticleDefinition;

class G4ionEffectiveCharge 
{

public:

  G4ionEffectiveCharge();

  virtual ~G4ionEffectiveCharge();

  G4double EffectiveChargeSquareRatio(
                           const G4ParticleDefinition* p,
                           const G4Material* material,
			         G4double kineticEnergy);

  G4double EffectiveCharge(const G4ParticleDefinition* p,
                           const G4Material* material,
			         G4double kineticEnergy);

private:

  // hide assignment operator
  G4ionEffectiveCharge & operator=(const G4ionEffectiveCharge &right);
  G4ionEffectiveCharge(const G4ionEffectiveCharge&);

  G4double                    chargeCorrection;
  G4double                    energyHighLimit;
  G4double                    energyLowLimit;
  G4double                    energyBohr;
  G4double                    massFactor;
  G4double                    minCharge;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ionEffectiveCharge::EffectiveChargeSquareRatio(
                           const G4ParticleDefinition* p,
                           const G4Material* material,
			         G4double kineticEnergy)
{
  G4double charge = EffectiveCharge(p,material,kineticEnergy)*chargeCorrection
                  / eplus;

  return charge*charge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
