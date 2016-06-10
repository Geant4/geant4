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
// $Id: G4RepleteEofM.cc $
//
//
//  This is the standard right-hand side for equation of motion.
//
//  08.04.2013 Peter Gumplinger
//
// -------------------------------------------------------------------

#include "G4RepleteEofM.hh"
#include "G4Field.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


G4RepleteEofM::G4RepleteEofM( G4Field* field, G4int nvar )
          : G4EquationOfMotion( field ), fNvar(nvar),
            fBfield(false), fEfield(false), fGfield(false), 
            fgradB(false), fSpin(false),
            charge(0.), mass(0.), magMoment(0.), spin(0.),
            ElectroMagCof(0.), omegac(0.), anomaly(0.),
            beta(0.), gamma(0.)
{
   fGfield = field->IsGravityActive();
}

G4RepleteEofM::~G4RepleteEofM()
{
}

void  
G4RepleteEofM::SetChargeMomentumMass(G4ChargeState particleCharge, // e+ units
                              G4double MomentumXc,
                              G4double particleMass)
{
   charge    = particleCharge.GetCharge();
   mass      = particleMass;
   magMoment = particleCharge.GetMagneticDipoleMoment();
   spin      = particleCharge.GetSpin();

   ElectroMagCof =  eplus*charge*c_light;
   omegac = (eplus/mass)*c_light;

   G4double muB = 0.5*eplus*hbar_Planck/(mass/c_squared);

   G4double g_BMT;
   if ( spin != 0. ) g_BMT = (std::abs(magMoment)/muB)/spin;
   else g_BMT = 2.;

   anomaly = (g_BMT - 2.)/2.;

   G4double E = std::sqrt(sqr(MomentumXc)+sqr(mass));
   beta  = MomentumXc/E;
   gamma = E/mass;
}

void
G4RepleteEofM::EvaluateRhsGivenB(const G4double y[],
                          const G4double Field[],
                          G4double dydx[] ) const
{

   // Components of y:
   //    0-2 dr/ds,
   //    3-5 dp/ds - momentum derivatives
   //    9-11 dSpin/ds = (1/beta) dSpin/dt - spin derivatives
   //
   // The BMT equation, following J.D.Jackson, Classical
   // Electrodynamics, Second Edition,
   // dS/dt = (e/mc) S \cross
   //              [ (g/2-1 +1/\gamma) B
   //               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta
   //               -(g/2-\gamma/(\gamma+1) \beta \cross E ]
   // where
   // S = \vec{s}, where S^2 = 1
   // B = \vec{B}
   // \beta = \vec{\beta} = \beta \vec{u} with u^2 = 1
   // E = \vec{E}
   //
   // Field[0,1,2] are the magnetic field components
   // Field[3,4,5] are the electric field components
   // Field[6,7,8] are the gravity  field components
   // The Field[] array may trivially be extended to 18 components
   // Field[ 9] == dB_x/dx; Field[10] == dB_y/dx; Field[11] == dB_z/dx
   // Field[12] == dB_x/dy; Field[13] == dB_y/dy; Field[14] == dB_z/dy
   // Field[15] == dB_x/dz; Field[16] == dB_y/dz; Field[17] == dB_z/dz

   G4double momentum_mag_square = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
   G4double inv_momentum_magnitude = 1.0 / std::sqrt( momentum_mag_square );

   G4double Energy = std::sqrt(momentum_mag_square + mass*mass);
   G4double inverse_velocity = Energy*inv_momentum_magnitude/c_light;

   G4double cof1 = ElectroMagCof*inv_momentum_magnitude;
   G4double cof2 = Energy/c_light;
   G4double cof3 = inv_momentum_magnitude*mass;

   dydx[0] = y[3]*inv_momentum_magnitude;       //  (d/ds)x = Vx/V
   dydx[1] = y[4]*inv_momentum_magnitude;       //  (d/ds)y = Vy/V
   dydx[2] = y[5]*inv_momentum_magnitude;       //  (d/ds)z = Vz/V

   dydx[3] = 0.;
   dydx[4] = 0.;
   dydx[5] = 0.;

   G4double field[18] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

   field[0] = Field[0];
   field[1] = Field[1];
   field[2] = Field[2];

   // Force due to B field - Field[0,1,2]

   if (fBfield) {
      if (charge != 0.) {
         dydx[3] += cof1*(y[4]*field[2] - y[5]*field[1]);
         dydx[4] += cof1*(y[5]*field[0] - y[3]*field[2]);
         dydx[5] += cof1*(y[3]*field[1] - y[4]*field[0]);
      }
   }

   // add force due to E field - Field[3,4,5]

   if (!fBfield) {
      field[3] = Field[0];
      field[4] = Field[1];
      field[5] = Field[2];
   } else {
      field[3] = Field[3];
      field[4] = Field[4];
      field[5] = Field[5];
   }

   if (fEfield) {
      if (charge != 0.) {
         dydx[3] += cof1*cof2*field[3];
         dydx[4] += cof1*cof2*field[4];
         dydx[5] += cof1*cof2*field[5];
      }
   }

   // add force due to gravity field - Field[6,7,8]

   if (!fBfield && !fEfield) {
      field[6] = Field[0];
      field[7] = Field[1];
      field[8] = Field[2];
   } else {
      field[6] = Field[6];
      field[7] = Field[7];
      field[8] = Field[8];
   }

   if (fGfield) {
      if (mass > 0.) {
         dydx[3] += field[6]*cof2*cof3/c_light;
         dydx[4] += field[7]*cof2*cof3/c_light;
         dydx[5] += field[8]*cof2*cof3/c_light;
      }
   }

   // add force due to ∇(µ⋅B) == (µ⋅∇)B when (∇xB) = 0

   if (!fBfield && !fEfield && !fGfield) {
      field[9]  = Field[0];
      field[10] = Field[1];
      field[11] = Field[2];
      field[12] = Field[3];
      field[13] = Field[4];
      field[14] = Field[5];
      field[15] = Field[6];
      field[16] = Field[7];
      field[17] = Field[8];
   } else {
      field[9]  = Field[9];
      field[10] = Field[10];
      field[11] = Field[11];
      field[12] = Field[12];
      field[13] = Field[13];
      field[14] = Field[14];
      field[15] = Field[15];
      field[16] = Field[16];
      field[17] = Field[17];
   }

   if (fgradB) {
      if (magMoment != 0.) {

         // field[ 9] == dB_x/dx; field[10] == dB_y/dx; field[11] == dB_z/dx
         // field[12] == dB_x/dy; field[13] == dB_y/dy; field[14] == dB_z/dy
         // field[15] == dB_x/dz; field[16] == dB_y/dz; field[17] == dB_z/dz

//         G4cout << "y[9]:  " << y[9] << " y[10]: " << y[10] << " y[11]: " << y[11] << G4endl;
//         G4cout << "field[9]:  " << field[9]  << " field[10]: " << field[10] << " field[11]: " << field[11] << G4endl;
//         G4cout << "field[12]: " << field[12] << " field[13]: " << field[13] << " field[14]: " << field[14] << G4endl;
//         G4cout << "field[15]: " << field[15] << " field[16]: " << field[16] << " field[17]: " << field[17] << G4endl;
//         G4cout << "inv_momentum_magnitdue: " << inv_momentum_magnitude << " Energy: " << Energy << G4endl;

         dydx[3] += magMoment*(y[9]*field[ 9]+y[10]*field[10]+y[11]*field[11])
                                                *inv_momentum_magnitude*Energy;
         dydx[4] += magMoment*(y[9]*field[12]+y[10]*field[13]+y[11]*field[14])
                                                *inv_momentum_magnitude*Energy;
         dydx[5] += magMoment*(y[9]*field[15]+y[10]*field[16]+y[11]*field[17])
                                                *inv_momentum_magnitude*Energy;

//         G4cout << "dydx[3,4,5] " << dydx[3] << " " << dydx[4] << " " << dydx[5] << G4endl;
      }
   }

   dydx[6] = 0.; //not used

   // Lab Time of flight
   dydx[7] = inverse_velocity;

   if (fNvar == 12) {
      dydx[ 8] = 0.; //not used

      dydx[ 9] = 0.;
      dydx[10] = 0.;
      dydx[11] = 0.;
   }

   if (fSpin) {
//      G4cout << "y[9,10,11]  " << y[9] << " " << y[10] << " " << y[11] << G4endl;
      G4ThreeVector BField(0.,0.,0.);
      if (fBfield) {
         G4ThreeVector F(field[0],field[1],field[2]);
         BField = F;
      }

      G4ThreeVector EField(0.,0.,0.);
      if (fEfield) {
         G4ThreeVector F(field[3],field[4],field[5]);
         EField = F;
      }

      EField /= c_light;

      G4ThreeVector u(y[3], y[4], y[5]);
      u *= inv_momentum_magnitude;

      G4double udb = anomaly*beta*gamma/(1.+gamma) * (BField * u);
      G4double ucb = (anomaly+1./gamma)/beta;
      G4double uce = anomaly + 1./(gamma+1.);

      G4ThreeVector Spin(y[9],y[10],y[11]);

      G4double pcharge;
      if (charge == 0.) pcharge = 1.;
      else pcharge = charge;

      G4ThreeVector dSpin(0.,0.,0);
      if (Spin.mag2() != 0.) {
         if (fBfield) {
           dSpin =
             pcharge*omegac*( ucb*(Spin.cross(BField))-udb*(Spin.cross(u)) );
         }
         if (fEfield) {
            dSpin -=
                                     // from Jackson
                                     // -uce*Spin.cross(u.cross(EField)) );
                                     // but this form has one less operation
             pcharge*omegac*( uce*(u*(Spin*EField) - EField*(Spin*u)) );
         }
      }

      dydx[ 9] = dSpin.x();
      dydx[10] = dSpin.y();
      dydx[11] = dSpin.z();

//      G4cout << "dydx[9,10,11] " << dydx[9] << " " << dydx[10] << " " << dydx[11] << G4endl;
   }

   return ;
}
