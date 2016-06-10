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
//      File name:     G4PionMinusField.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------

#include "G4PionMinusField.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NucleiProperties.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4PionMinus.hh"
#include "G4HadTmpUtil.hh"
#include "G4Pow.hh"

G4PionMinusField::G4PionMinusField(G4V3DNucleus * nucleus, G4double coeff)
  : G4VNuclearField(nucleus)
{
  theCoeff = coeff;
}


G4PionMinusField::~G4PionMinusField()
{ }

G4double G4PionMinusField::GetField(const G4ThreeVector & aPosition)
{
// Field is 0 out of the nucleus!
  if(aPosition.mag() >= radius) return 0.0;

  G4int A = theNucleus->GetMassNumber();
  G4int Z = theNucleus->GetCharge();
  G4double pionMinusMass = G4PionMinus::PionMinus()->GetPDGMass();
  G4double bindingEnergy = G4NucleiProperties::GetBindingEnergy(G4lrint(A), G4lrint(Z));
  G4double nucleusMass = Z*proton_mass_c2+(A-Z)*neutron_mass_c2+bindingEnergy;
  G4double reducedMass = pionMinusMass*nucleusMass/(pionMinusMass+nucleusMass);

  G4double density = A*theNucleus->GetNuclearDensity()->GetDensity(aPosition);
  G4double nucleonMass = (proton_mass_c2+neutron_mass_c2)/2;

  return 2.*pi*hbarc*hbarc/reducedMass*(1+pionMinusMass/nucleonMass)*theCoeff*density + GetBarrier();
}


G4double G4PionMinusField::GetBarrier()
{
  G4int A = theNucleus->GetMassNumber();
  G4int Z = theNucleus->GetCharge();
  G4double coulombBarrier = (1.44/1.14) * MeV * Z / (1.0 + G4Pow::GetInstance()->Z13(A));
  return -coulombBarrier;
}
















