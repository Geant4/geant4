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
// $Id: G4EqEMFieldWithEDM.cc 95822 2016-02-26 08:04:51Z gcosmo $
//
//
//  This is the standard right-hand side for equation of motion.
//
//  19.02.2009 Kevin Lynch, based on G4EqEMFieldWithSpin
//  06.11.2009 Hiromi Iinuma see:
//  http://hypernews.slac.stanford.edu/HyperNews/geant4/get/emfields/161.html
//
// -------------------------------------------------------------------

#include "G4EqEMFieldWithEDM.hh"
#include "G4ElectroMagneticField.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4EqEMFieldWithEDM::G4EqEMFieldWithEDM(G4ElectroMagneticField *emField )
  : G4EquationOfMotion( emField ), charge(0.), mass(0.), magMoment(0.),
    spin(0.), fElectroMagCof(0.), fMassCof(0.), omegac(0.),
    anomaly(0.0011659208), eta(0.), beta(0.), gamma(0.)
{
}

G4EqEMFieldWithEDM::~G4EqEMFieldWithEDM()
{
} 

void  
G4EqEMFieldWithEDM::SetChargeMomentumMass(G4ChargeState particleCharge,
                                          G4double MomentumXc,
                                          G4double particleMass)
{
   charge   = particleCharge.GetCharge();
   mass      = particleMass;
   magMoment = particleCharge.GetMagneticDipoleMoment();
   spin      = particleCharge.GetSpin();

   fElectroMagCof =  eplus*charge*c_light;
   fMassCof = mass*mass;

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
G4EqEMFieldWithEDM::EvaluateRhsGivenB(const G4double y[],
                                       const G4double Field[],
                                             G4double dydx[] ) const
{

   // Components of y:
   //    0-2 dr/ds,
   //    3-5 dp/ds - momentum derivatives
   //    9-11 dSpin/ds = (1/beta) dSpin/dt - spin derivatives

   // The BMT equation, following J.D.Jackson, Classical
   // Electrodynamics, Second Edition, with additions for EDM
   // evolution from 
   // M.Nowakowski, et.al. Eur.J.Phys.26, pp 545-560, (2005)
   // or
   // Silenko, Phys.Rev.ST Accel.Beams 9:034003, (2006)

   // dS/dt = (e/m) S \cross 
   // MDM:         [ (g/2-1 +1/\gamma) B
   //               -(g/2-1)\gamma/(\gamma+1) (\beta \cdot B)\beta 
   //               -(g/2-\gamma/(\gamma+1) \beta \cross E 
   //
   // EDM:        eta/2( E - gamma/(gamma+1) \beta (\beta \cdot E)
   //                    + \beta \cross B ) ]
   //
   // where
   // S = \vec{s}, where S^2 = 1
   // B = \vec{B}
   // \beta = \vec{\beta} = \beta \vec{u} with u^2 = 1
   // E = \vec{E}

   G4double pSquared = y[3]*y[3] + y[4]*y[4] + y[5]*y[5] ;

   G4double Energy   = std::sqrt( pSquared + fMassCof );
   G4double cof2     = Energy/c_light ;

   G4double pModuleInverse  = 1.0/std::sqrt(pSquared) ;

   G4double inverse_velocity = Energy * pModuleInverse / c_light;

   G4double cof1     = fElectroMagCof*pModuleInverse ;

   dydx[0] = y[3]*pModuleInverse ;                         
   dydx[1] = y[4]*pModuleInverse ;                         
   dydx[2] = y[5]*pModuleInverse ;                        

   dydx[3] = cof1*(cof2*Field[3] + (y[4]*Field[2] - y[5]*Field[1])) ;
   
   dydx[4] = cof1*(cof2*Field[4] + (y[5]*Field[0] - y[3]*Field[2])) ; 
 
   dydx[5] = cof1*(cof2*Field[5] + (y[3]*Field[1] - y[4]*Field[0])) ;  
   
   dydx[6] = dydx[8] = 0.;//not used

   // Lab Time of flight
   dydx[7] = inverse_velocity;
   
   G4ThreeVector BField(Field[0],Field[1],Field[2]);
   G4ThreeVector EField(Field[3],Field[4],Field[5]);

   EField /= c_light;

   G4ThreeVector u(y[3], y[4], y[5]);
   u *= pModuleInverse;

   G4double udb = anomaly*beta*gamma/(1.+gamma) * (BField * u);
   G4double ucb = (anomaly+1./gamma)/beta;
   G4double uce = anomaly + 1./(gamma+1.);
   G4double ude = beta*gamma/(1.+gamma)*(EField*u);

   G4ThreeVector Spin(y[9],y[10],y[11]);

   G4double pcharge;
   if (charge == 0.) pcharge = 1.;
   else pcharge = charge;

   G4ThreeVector dSpin(0.,0.,0.);
   if (Spin.mag2() != 0.) {
      dSpin =
         pcharge*omegac*( ucb*(Spin.cross(BField))-udb*(Spin.cross(u))
                                 // from Jackson
                                 // -uce*Spin.cross(u.cross(EField)) )
                                 // but this form has one less operation
                         - uce*(u*(Spin*EField) - EField*(Spin*u))
                         + eta/2.*(Spin.cross(EField) - ude*(Spin.cross(u))
                                 // +Spin.cross(u.cross(Bfield))
                         + (u*(Spin*BField) - BField*(Spin*u)) ) );
   }
      
   dydx[ 9] = dSpin.x();
   dydx[10] = dSpin.y();
   dydx[11] = dSpin.z();

   return ;
}
