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
// $Id: G4ionEffectiveCharge.hh 100363 2016-10-19 09:24:47Z gcosmo $
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
#include "G4ParticleDefinition.hh"

class G4Material;
class G4Pow;

class G4ionEffectiveCharge 
{

public:

  explicit G4ionEffectiveCharge();

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
  G4ionEffectiveCharge & operator=(const G4ionEffectiveCharge &right) = delete;
  G4ionEffectiveCharge(const G4ionEffectiveCharge&) = delete;

  static const G4double       inveplus;
  G4Pow*                      g4calc;

  const G4ParticleDefinition* lastPart;
  const G4Material*           lastMat;
  G4double                    lastKinEnergy;

  G4double                    chargeCorrection;
  G4double                    effCharge;

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
  G4double charge = effCharge;
  if( kineticEnergy != lastKinEnergy || material != lastMat || p != lastPart) {
    charge = EffectiveCharge(p,material,kineticEnergy);
  }
  charge *= chargeCorrection*inveplus;

  return charge*charge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
