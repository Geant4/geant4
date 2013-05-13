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
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4KM_OpticalEqRhs.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#include "G4KM_OpticalEqRhs.hh"
#include "G4PhysicalConstants.hh"
#include "G4NucleiProperties.hh"
#include "G4VNuclearDensity.hh"
#include "G4HadTmpUtil.hh"

G4KM_OpticalEqRhs::G4KM_OpticalEqRhs(G4KM_DummyField *field,
				     G4V3DNucleus * nucleus) :
  G4Mag_EqRhs(field), theNucleus(nucleus)
{
  theFactor = 0;
  theMass = 0;
}


void G4KM_OpticalEqRhs::SetFactor(G4double mass, G4double opticalParameter)
{
  G4double A = theNucleus->GetMassNumber();
  G4double Z = theNucleus->GetCharge();
  G4double bindingEnergy = G4NucleiProperties::GetBindingEnergy(G4lrint(A), G4lrint(Z));
  G4double nucleusMass = Z*proton_mass_c2+(A-Z)*neutron_mass_c2+bindingEnergy;
  G4double reducedMass = mass*nucleusMass/(mass+nucleusMass);

  G4double nucleonMass = (proton_mass_c2+neutron_mass_c2)/2;

// _factor in (MeV*fermi)*fermi/MeV = fermi*fermi  -- need to have A as density normalized to 1
  theFactor = 2*pi*hbarc*hbarc*(1+mass/nucleonMass)* opticalParameter/reducedMass * A;

  theMass = mass;
}


void G4KM_OpticalEqRhs::EvaluateRhsGivenB(const G4double y[], const G4double *,
					  G4double dydx[]) const
{
  G4double yMod = std::sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
  G4double e = std::sqrt(theMass*theMass+y[3]*y[3]+y[4]*y[4]+y[5]*y[5]);
  dydx[0] = c_light*y[3]/e;   //
  dydx[1] = c_light*y[4]/e;   //  dq/dt=dH/dp = c*p/e
  dydx[2] = c_light*y[5]/e;   // 

// V=K*rho(r) ==> dydx[3] = -dV/dr*dr/dx = -K*d(rho)/dr*dr/dx.
// Idem for dydx[4] and dydx[5]

  const G4VNuclearDensity * nuclearDensity=theNucleus->GetNuclearDensity();

  G4ThreeVector pos(y[0],y[1],y[2]);
  G4double deriv = theFactor*nuclearDensity->GetDeriv(pos);

  dydx[3] = yMod == 0 ? 0 : -deriv*y[0]/yMod*c_light;
  dydx[4] = yMod == 0 ? 0 : -deriv*y[1]/yMod*c_light;
  dydx[5] = yMod == 0 ? 0 : -deriv*y[2]/yMod*c_light;
}

// Here by design, but it is unnecessary for nuclear fields
void G4KM_OpticalEqRhs::SetChargeMomentumMass(G4ChargeState,G4double ,G4double )
{ 
}
