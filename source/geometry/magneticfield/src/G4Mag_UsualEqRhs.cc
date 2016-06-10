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
// $Id: G4Mag_UsualEqRhs.cc 69699 2013-05-13 08:50:30Z gcosmo $
//
//
//  This is the 'standard' right-hand side for the equation of motion
//  of a charged particle in a magnetic field.
//
//  Initial version: J. Apostolakis, January 13th, 1997
//
// --------------------------------------------------------------------

#include "G4Mag_UsualEqRhs.hh"
#include "G4MagneticField.hh"

#include "globals.hh"    // For DBL_MAX

G4Mag_UsualEqRhs::G4Mag_UsualEqRhs( G4MagneticField* MagField )
  : G4Mag_EqRhs( MagField ) {}

G4Mag_UsualEqRhs::~G4Mag_UsualEqRhs() {}

void
G4Mag_UsualEqRhs::EvaluateRhsGivenB( const G4double y[],
			             const G4double B[3],
				           G4double dydx[] ) const
{
   G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
   G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );

   G4double cof = FCof()*inv_momentum_magnitude;

   dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

   return ;
}

void
G4Mag_UsualEqRhs::
 SetChargeMomentumMass( G4ChargeState particleCharge,
                        G4double MomentumXc,
			G4double mass)

{
   G4Mag_EqRhs::SetChargeMomentumMass( particleCharge, MomentumXc, mass);
}
