// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Mag_UsualEqRhs.cc,v 1.1 1999-01-07 16:07:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  This is the standard right-hand side for equation of motion.
//
//    The only case another is required is when using a moving reference
//     frame ... or extending the class to include additional Forces,
//     eg an electric field
//
//            J. Apostolakis, January 13th, 1997
//
#include "G4Mag_UsualEqRhs.hh"

void
G4Mag_UsualEqRhs::EvaluateRhsGivenB( const G4double y[],
			      const G4double B[3],
				    G4double dydx[] ) const
{
   G4double velocity_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
   G4double inv_velocity_magnitude = 1.0 / sqrt( velocity_mag_square );

   dydx[0] = y[3] * inv_velocity_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4] * inv_velocity_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5] * inv_velocity_magnitude;       //  (d/ds)z = Vz/V
   dydx[3] = FCof()*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = FCof()*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = FCof()*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

   return ;
}
