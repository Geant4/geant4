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
//  This is the right-hand side for equation of motion for a
//  massive particle in a gravitational field.
//
// History:
// - 14.06.11 P.Gumplinger, Created.
// -------------------------------------------------------------------
// Adopted from G4EqMagElectricField.hh
//
// Thanks to Peter Fierlinger (PSI) and
// A. Capra and A. Fontana (INFN Pavia)
// -------------------------------------------------------------------

#include "G4EqGravityField.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"

void
G4EqGravityField::SetChargeMomentumMass(G4ChargeState,
                                        G4double,
                                        G4double particleMass )
{
  fMass = particleMass;
}

void
G4EqGravityField::EvaluateRhsGivenB(const G4double y[],
                                    const G4double G[],
                                    G4double dydx[] ) const
{

  // Components of y:
  //    0-2 dr/ds,
  //    3-5 dp/ds - momentum derivatives

  G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
  G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );

  G4double Energy = std::sqrt(momentum_mag_square + fMass*fMass);
  G4double cof2 = Energy/c_light;
  G4double cof1 = inv_momentum_magnitude*fMass;
  G4double inverse_velocity = Energy*inv_momentum_magnitude/c_light;

  dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
  dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
  dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

  dydx[3] = G[0]*cof1*cof2/c_light;
  dydx[4] = G[1]*cof1*cof2/c_light;            //  m*g
  dydx[5] = G[2]*cof1*cof2/c_light;

  // Lab Time of flight
  dydx[7] = inverse_velocity;

  return;
}
