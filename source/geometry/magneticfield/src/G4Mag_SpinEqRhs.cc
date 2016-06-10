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
// $Id: G4Mag_SpinEqRhs.cc 69970 2013-05-21 10:14:30Z gcosmo $
//
// This is the standard right-hand side for equation of motion.
// This version of the right-hand side includes the three components
// of the particle's spin.
//
//            J. Apostolakis, February 8th, 1999
//            P. Gumplinger,  February 8th, 1999
//            D. Cote-Ahern, P. Gumplinger,  April 11th, 2001
//
// --------------------------------------------------------------------

#include "G4Mag_SpinEqRhs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

G4Mag_SpinEqRhs::G4Mag_SpinEqRhs( G4MagneticField* MagField )
  : G4Mag_EqRhs( MagField ), omegac(0.), anomaly(0.0011659208),
    pcharge(0.), E(0.), gamma(0.), beta(0.)
{
}

G4Mag_SpinEqRhs::~G4Mag_SpinEqRhs()
{
}

void
G4Mag_SpinEqRhs::SetChargeMomentumMass(G4double particleCharge, // in e+ units
                                       G4double MomentumXc,
                                       G4double particleMass)
{
   //  To set fCof_val 
   G4Mag_EqRhs::SetChargeMomentumMass(particleCharge, MomentumXc, particleMass);

   omegac = (eplus/particleMass)*c_light;

   pcharge = particleCharge;

   E = std::sqrt(sqr(MomentumXc)+sqr(particleMass));
   beta  = MomentumXc/E;
   gamma = E/particleMass;

   G4double neutronAnomaly = -2.9156797;
   if (pcharge==0.) SetAnomaly(neutronAnomaly);
}

void
G4Mag_SpinEqRhs::EvaluateRhsGivenB( const G4double y[],
                                    const G4double B[3],
                                          G4double dydx[] ) const
{
   G4double momentum_mag_square = sqr(y[3]) + sqr(y[4]) + sqr(y[5]);
   G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );
   G4double cof = FCof()*inv_momentum_magnitude;

   dydx[0] = y[3] * inv_momentum_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4] * inv_momentum_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5] * inv_momentum_magnitude;       //  (d/ds)z = Vz/V

   if (pcharge == 0.) {
      dydx[3] = 0.;
      dydx[4] = 0.;
      dydx[5] = 0.;
   } else {
      dydx[3] = cof*(y[4]*B[2] - y[5]*B[1]) ;   // Ax = a*(Vy*Bz - Vz*By)
      dydx[4] = cof*(y[5]*B[0] - y[3]*B[2]) ;   // Ay = a*(Vz*Bx - Vx*Bz)
      dydx[5] = cof*(y[3]*B[1] - y[4]*B[0]) ;   // Az = a*(Vx*By - Vy*Bx)
   }

   G4ThreeVector u(y[3], y[4], y[5]);
   u *= inv_momentum_magnitude; 

   G4ThreeVector BField(B[0],B[1],B[2]);

   G4double udb = anomaly*beta*gamma/(1.+gamma) * (BField * u); 
   G4double ucb = (anomaly+1./gamma)/beta;

   // Initialise the values of dydx that we do not update.
   dydx[6] = dydx[7] = dydx[8] = 0.0;

   G4ThreeVector Spin(y[9],y[10],y[11]);

   G4ThreeVector dSpin;

   if (pcharge == 0.) {
      // dSpin = (3.8260837/2.)*omegac*(Spin.cross(BField));
      dSpin = omegac*(ucb*(Spin.cross(BField))-udb*(Spin.cross(u)));
   } else {
      dSpin = pcharge*omegac*(ucb*(Spin.cross(BField))-udb*(Spin.cross(u)));
   }

   dydx[ 9] = dSpin.x();
   dydx[10] = dSpin.y();
   dydx[11] = dSpin.z();

   return ;
}
