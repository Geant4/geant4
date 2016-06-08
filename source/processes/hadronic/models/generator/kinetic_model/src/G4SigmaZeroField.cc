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
//      File name:     G4SigmaZeroField.cc
//
//      Author:        Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
// 
//      Creation date: 5 June 2000
// -------------------------------------------------------------------
#include "G4SigmaZeroField.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4SigmaZero.hh"

G4SigmaZeroField::G4SigmaZeroField(G4V3DNucleus * nucleus, G4double coeff)
  : G4VNuclearField(nucleus)
{
  theCoeff = coeff;
}


G4SigmaZeroField::~G4SigmaZeroField()
{ }


const G4SigmaZeroField & G4SigmaZeroField::operator=(const G4SigmaZeroField & right)
{
  G4Exception("G4SigmaZeroField::operator= meant not to be accessible");
  return *this;
}


G4int G4SigmaZeroField::operator==(const G4SigmaZeroField & right) const
{
  G4Exception("G4SigmaZeroField::operator== meant not to be accessible");
  return 0;
}


G4int G4SigmaZeroField::operator!=(const G4SigmaZeroField & right) const
{
  G4Exception("G4SigmaZeroField::operator!= meant not to be accessible");
  return 1;
}



G4double G4SigmaZeroField::GetField(const G4ThreeVector & aPosition)
{
// Field is 0 out of the nucleus!
  if(aPosition.mag() >= radius) return 0.0;

  G4double sigmaZeroMass = G4SigmaZero::SigmaZero()->GetPDGMass();

  G4double A = theNucleus->GetMassNumber();
  G4double Z = theNucleus->GetCharge();
  G4double bindingEnergy = G4NucleiPropertiesTable::GetBindingEnergy(Z, A);
  G4double nucleusMass = Z*proton_mass_c2+(A-Z)*neutron_mass_c2+bindingEnergy;
  G4double reducedMass = sigmaZeroMass*nucleusMass/(sigmaZeroMass+nucleusMass);

  const G4VNuclearDensity * nuclearDensity=theNucleus->GetNuclearDensity();
  G4double density = nuclearDensity->GetDensity(aPosition);

  return -2.*pi*hbarc*hbarc/reducedMass*(2.0)*theCoeff*density;
}


G4double G4SigmaZeroField::GetBarrier()
{
  return 0.;
}

