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
// $Id: G4Mag_UsualEqRhs.cc,v 1.6 2002-05-03 16:01:52 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  This is the 'standard' right-hand side for the equation of motion
//    of a charged particle in a magnetic field.
//
//  Initial version: J. Apostolakis, January 13th, 1997
//  Modified:  
//             J. Apostolakis, April 4th, 2002: for speedup
//
#include "G4Mag_UsualEqRhs.hh"

void
G4Mag_UsualEqRhs::EvaluateRhsGivenB( const G4double y[],
			      const G4double B[3],
				    G4double dydx[] ) const
{
   G4double cof; 
   register double crossX, crossY, crossZ; 
   G4double inv_momentum_magnitude;

   // Original was
   //    inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );

   G4double momentum_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
   G4double epsil_mom_sq= momentum_mag_square * sqr(fInvCurrentMomentumXc) - 1.0;
   // G4double inv_momentum_mag_1stOrder= fInvCurrentMomentumXc * 
     // ( 1.0 - 0.5 * epsil_mom_sq );

   G4double inv_momentum_mag_2ndOrder= fInvCurrentMomentumXc * 
     ( 1.0 - 0.5 * epsil_mom_sq + 0.375 * epsil_mom_sq * epsil_mom_sq );

   inv_momentum_magnitude= inv_momentum_mag_2ndOrder;
   //******************************************************
#if 0   
   G4double currentMomentum= 1.0 / fInvCurrentMomentumXc;
   G4double inv_momentum_mag_other= currentMomentum / momentum_mag_square;

   G4double inv_momentum_mag_other_1stOrder= inv_momentum_mag_other *
     ( 1.0 + 0.5 * epsil_mom_sq );

   inv_momentum_magnitude= inv_momentum_mag_other_1stOrder;
   //******************************************************
#endif 

   cof = FCof()*inv_momentum_magnitude;

   crossX= y[4]*B[2] - y[5]*B[1] ;   //  Cx = Vy*Bz - Vz*By
   crossY= y[5]*B[0] - y[3]*B[2] ;   // Cy = Vz*Bx - Vx*Bz
   crossZ= y[3]*B[1] - y[4]*B[0] ;   // Cz = Vx*By - Vy*Bx


   dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

   dydx[3] = cof*crossX ; // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*crossY ; // (y[5]*B[0] - y[3]*B[2]); ie Ay= a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*crossZ ; // (y[3]*B[1] - y[4]*B[0]); ie Az= a*(Vx*By - Vy*Bx)

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
