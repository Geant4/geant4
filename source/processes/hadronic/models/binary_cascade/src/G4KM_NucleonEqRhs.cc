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
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4KM_NucleonEqRhs.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#include "G4KM_NucleonEqRhs.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4VNuclearDensity.hh"


G4KM_NucleonEqRhs::G4KM_NucleonEqRhs(G4KM_DummyField *field,
				     G4V3DNucleus * nucleus) :
  G4Mag_EqRhs(field), theNucleus(nucleus)
{
  theMass = 0.;
  A = theNucleus->GetMassNumber();
  factor = hbarc*hbarc*pow(3.*pi2*A,2./3.)/3.;
}


void G4KM_NucleonEqRhs::EvaluateRhsGivenB(const G4double y[],
					  const G4double *,
					  G4double dydx[]) const
{
  G4double yMod = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
  G4double e = sqrt(theMass*theMass+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);

// y[0..2] is position
// y[3..5] is momentum (and not mom.direction)
    
  dydx[0] = c_light*y[3]/e;   //
  dydx[1] = c_light*y[4]/e;   //  dq/dt=dH/dp = c*p/e
  dydx[2] = c_light*y[5]/e;   // 

/*
 * // debug
 *  G4cout << " Nucleon RHS : 0..2(dpos/dt) " << 
 *       dydx[0] << " " <<
 *       dydx[1] << " " <<
 *       dydx[2] << " " << G4endl;
 */

      
// V=K*rho(r) ==> dydx[3] = -dV/dr*dr/dx = -K*d(rho)/dr*dr/dx.
// GF should be V=K*rho(r) ==> dydx[3] = -dV/dr*dr/dx = -K*d(rho)/dr*dr/dt
// GF  and dV/dt = dE/dt ==> dp/dt = dE/dt * dp/dE = dE/dt *e/p
// Idem for dydx[4] and dydx[5]

  G4ThreeVector pos(y[0],y[1],y[2]);

  const G4VNuclearDensity * nuclearDensity=theNucleus->GetNuclearDensity();

// do not check for theMass != 0 : it is an error and core dump will signal it

  G4double density=  nuclearDensity->GetDensity(pos);
  G4double deriv(0);
  if (density > 0 ) deriv = (factor/theMass)*
			pow(density, -1./3.)*nuclearDensity->GetDeriv(pos);

//  dydx[3] = yMod == 0 ? 0 : -deriv*y[0]/yMod;
//  dydx[4] = yMod == 0 ? 0 : -deriv*y[1]/yMod;
//  dydx[5] = yMod == 0 ? 0 : -deriv*y[2]/yMod;
  dydx[3] = yMod == 0 ? 0 : deriv*y[0]/yMod*c_light;
  dydx[4] = yMod == 0 ? 0 : deriv*y[1]/yMod*c_light;
  dydx[5] = yMod == 0 ? 0 : deriv*y[2]/yMod*c_light;


/*
 * // debug
 * G4cout << " Nucleon RHS : 3..5(dE/dt) " << 
 *       dydx[3] << " " <<
 *       dydx[4] << " " <<
 *       dydx[5] << " " << G4endl;
 */
}




