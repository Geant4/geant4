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
//      File name:     G4PionPlusField.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#include "G4PionPlusField.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4PionPlus.hh"
#include "G4HadTmpUtil.hh"
#include "G4Pow.hh"

G4PionPlusField::G4PionPlusField(G4V3DNucleus * nucleus, G4double coeff)
  : G4VNuclearField(nucleus)
{
  theCoeff = coeff;
}


G4PionPlusField::~G4PionPlusField()
{ }

G4double G4PionPlusField::GetField(const G4ThreeVector & aPosition)
{
// Field is 0 out of the nucleus!
  if(aPosition.mag() >= radius) return 0.0;

  G4double pionPlusMass = G4PionPlus::PionPlus()->GetPDGMass();

  G4int A = theNucleus->GetMassNumber();
  G4int Z = theNucleus->GetCharge();
  G4double bindingEnergy = G4NucleiProperties::GetBindingEnergy(A, Z);
  G4double nucleusMass = Z*proton_mass_c2+(A-Z)*neutron_mass_c2+bindingEnergy;
  G4double reducedMass = pionPlusMass*nucleusMass/(pionPlusMass+nucleusMass);

  G4double density = A*theNucleus->GetNuclearDensity()->GetDensity(aPosition);
  G4double nucleonMass = (proton_mass_c2+neutron_mass_c2)/2;

  return 2.*pi*hbarc*hbarc/reducedMass*(1+pionPlusMass/nucleonMass)*theCoeff*density + GetBarrier();
}


G4double G4PionPlusField::GetBarrier()
{
  G4int A = theNucleus->GetMassNumber();
  G4int Z = theNucleus->GetCharge();
  G4double coulombBarrier = (1.44/1.14) * MeV * Z / (1.0 + G4Pow::GetInstance()->Z13(A));
  return coulombBarrier;
}
