//
//
//  This is the standard right-hand side for equation of motion.
//        This version of the right-hand side includes
//        the three components of the particle's spin.
//
//            J. Apostolakis, February 8th, 1999
//            P. Gumplinger,  February 8th, 1999
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

   // for testing only
   anomaly = 0.0; 
}

void
G4Mag_SpinEqRhs::EvaluateRhsGivenB( const G4double y[],
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

   G4double beta_squared = velocity_mag_square/c_squared;
   G4double beta  = sqrt(beta_squared);
  
   G4double gamma;

   if (beta < 1.0){
      gamma = 1. / sqrt( 1. - beta_squared);
   } else {
      beta = 1.0; 
      gamma = DBL_MAX; 
   }

   G4ThreeVector u;
   u.setX(inv_velocity_magnitude*y[3]);
   u.setY(inv_velocity_magnitude*y[4]);
   u.setZ(inv_velocity_magnitude*y[5]);

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
