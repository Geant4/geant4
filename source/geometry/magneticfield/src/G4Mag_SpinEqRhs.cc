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
//
//  This is the standard right-hand side for equation of motion.
//        This version of the right-hand side includes
//        the three components of the particle's spin.
//
//            J. Apostolakis, February 8th, 1999
//            P. Gumplinger,  February 8th, 1999
//            D. Cote-Ahern, P. Gumplinger,  April 11th, 2001
//
#include "G4Mag_SpinEqRhs.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

void
G4Mag_SpinEqRhs::SetChargeMomentumMass(G4double particleCharge, // in e+ units
                                       G4double MomentumXc,
                                       G4double mass)
{
   //  To set fCof_val 
   G4Mag_EqRhs::SetChargeMomentumMass(particleCharge, MomentumXc, mass);

   omegac = 0.105658387*GeV/mass * 2.837374841e-3*(rad/cm/kilogauss);
   anomaly = 1.165923e-3;
   ParticleCharge = particleCharge;

   E = sqrt(sqr(MomentumXc)+sqr(mass));
   beta  = MomentumXc/E;
   gamma = E/mass;

}

void
G4Mag_SpinEqRhs::EvaluateRhsGivenB( const G4double y[],
			            const G4double B[3],
				    G4double dydx[] ) const
{
   G4double momentum_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
   G4double inv_momentum_magnitude = 1.0 / sqrt( momentum_mag_square );
   G4double cof = FCof()*inv_momentum_magnitude;

   dydx[0] = y[3] * inv_momentum_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4] * inv_momentum_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5] * inv_momentum_magnitude;       //  (d/ds)z = Vz/V
   dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
   dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
   dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)

   G4ThreeVector u(y[3], y[4], y[5]);
   u *= inv_momentum_magnitude; 

   G4ThreeVector BField(B[0],B[1],B[2]);

   G4double udb = anomaly*beta*gamma/(1.+gamma) * (BField * u); 
   G4double ucb = (anomaly+1./gamma)/beta;

   G4ThreeVector Spin(y[9],y[10],y[11]);

   G4ThreeVector dSpin;

   dSpin = ParticleCharge*omegac*(ucb*(Spin.cross(BField))-udb*(Spin.cross(u)));

   dydx[ 9] = dSpin.x();
   dydx[10] = dSpin.y();
   dydx[11] = dSpin.z();

   return ;
}
