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
//      File name:     G4PionZeroField.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------
#include "G4PionZeroField.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4PionZero.hh"


G4PionZeroField::G4PionZeroField(G4V3DNucleus * nucleus, G4double coeff)
  : G4VNuclearField(nucleus)
{
  theCoeff = coeff;
}


G4PionZeroField::~G4PionZeroField()
{ }


const G4PionZeroField & G4PionZeroField::operator=(const G4PionZeroField & right)
{
  G4Exception("G4PionZeroField::operator= meant not to be accessible");
  return *this;
}


G4int G4PionZeroField::operator==(const G4PionZeroField & right) const
{
  G4Exception("G4PionZeroField::operator== meant not to be accessible");
  return 0;
}


G4int G4PionZeroField::operator!=(const G4PionZeroField & right) const
{
  G4Exception("G4PionZeroField::operator!= meant not to be accessible");
  return 1;
}


G4double G4PionZeroField::GetField(const G4ThreeVector & aPosition)
{
// Field is 0 out of the nucleus!
  if(aPosition.mag() >= radius) return 0.0;

  G4double pionZeroMass = G4PionZero::PionZero()->GetPDGMass();
  G4double A = theNucleus->GetMassNumber();
  G4double Z = theNucleus->GetCharge();

  G4double bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(Z, A);
  G4double nucleusMass = Z*proton_mass_c2+(A-Z)*neutron_mass_c2+bindingEnergy;
  G4double reducedMass = pionZeroMass*nucleusMass/(pionZeroMass+nucleusMass);


  G4double density = A*theNucleus->GetNuclearDensity()->GetDensity(aPosition);
  G4double nucleonMass = (proton_mass_c2+neutron_mass_c2)/2;

  return 2.*pi*hbarc*hbarc/reducedMass*(1+pionZeroMass/nucleonMass)*theCoeff*density;
}

G4double G4PionZeroField::GetBarrier()
{
  return 0;
}
