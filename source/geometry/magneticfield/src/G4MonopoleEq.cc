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
// $Id: G4MonopoleEq.cc 69699 2013-05-13 08:50:30Z gcosmo $
//
//
//  This is the right-hand side for equation of motion for a 
//   magnetic charge in a combined Electro-Magnetic field
//
//  d(p_c)/ds=g{c-energyB_ - p_c x E}/pc
//
//  17.11.09   V.Grichine
//
// -------------------------------------------------------------------

#include "G4MonopoleEq.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

void  
G4MonopoleEq::SetChargeMomentumMass(G4ChargeState particleCharge, // e+ units
		                    G4double,
                                    G4double particleMass)
{
  G4double pcharge = particleCharge.GetCharge();
  fElectroMagCof =  eplus*pcharge;  // no *c_light as for ususal q
  fElectroMagCof /= 2*fine_structure_const;

  fMassCof = particleMass*particleMass ; 
}



void
G4MonopoleEq::EvaluateRhsGivenB(const G4double y[],
			                const G4double Field[],
				              G4double dydx[] ) const
{

   // Components of y:
   //    0-2 dr/ds, 
   //    3-5 d(pc)/ds - momentum derivatives 

   G4double pSquared = y[3]*y[3] + y[4]*y[4] + y[5]*y[5] ;

   G4double Energy   = std::sqrt( pSquared + fMassCof );
   G4double cof2     = Energy*c_light ;

   G4double pModuleInverse  = 1.0/std::sqrt(pSquared) ;

   //  G4double inverse_velocity = Energy * c_light * pModuleInverse;
   G4double inverse_velocity = Energy * pModuleInverse / c_light;

   G4double cof1     = fElectroMagCof*pModuleInverse ;

   //  G4double vDotE = y[3]*Field[3] + y[4]*Field[4] + y[5]*Field[5] ;

   dydx[0] = y[3]*pModuleInverse ;                         
   dydx[1] = y[4]*pModuleInverse ;                         
   dydx[2] = y[5]*pModuleInverse ;                        

   dydx[3] = cof1*(cof2*Field[0] - (y[4]*Field[5] - y[5]*Field[4])) ;
   
   dydx[4] = cof1*(cof2*Field[1] - (y[5]*Field[3] - y[3]*Field[5])) ; 
 
   dydx[5] = cof1*(cof2*Field[2] - (y[3]*Field[4] - y[4]*Field[3])) ;  

   dydx[6] = 0.;//not used

   // Lab Time of flight
   dydx[7] = inverse_velocity;
   return ;
}
