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
// $Id: G4Mag_UsualEqRhs.cc,v 1.8 2002-11-09 02:58:14 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  This is the 'standard' right-hand side for the equation of motion
//    of a charged particle in a magnetic field.
//
//  Initial version: J. Apostolakis, January 13th, 1997
//  Modified:        J. Apostolakis, November 9th, 2002
//
#include "G4Mag_UsualEqRhs.hh"

void
G4Mag_UsualEqRhs::EvaluateRhsGivenB( const G4double y[],
			      const G4double B[3],
				    G4double dydx[] ) const
{
   G4double momentum_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
   G4double inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );

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
 SetChargeMomentumMass( G4double particleCharge, // in e+ units
			                 G4double MomentumXc,
			                 G4double mass)

{
   fInvCurrentMomentumXc= 1.0 / MomentumXc;
   G4Mag_EqRhs::SetChargeMomentumMass( particleCharge, MomentumXc, mass);
}
