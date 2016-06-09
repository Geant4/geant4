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
// $Id: G4Be8FermiFragment.cc,v 1.5 2003/12/02 12:22:03 larazb Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4Be8FermiFragment.hh"
#include "G4IonTable.hh"
#include "G4HadronicException.hh"

G4Be8FermiFragment::G4Be8FermiFragment()
{
}

G4Be8FermiFragment::G4Be8FermiFragment(const G4Be8FermiFragment &) : G4UnstableFermiFragment()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4Be8FermiFragment::copy_constructor meant to not be accessable");
}


G4Be8FermiFragment::~G4Be8FermiFragment()
{
}


const G4Be8FermiFragment & G4Be8FermiFragment::operator=(const G4Be8FermiFragment &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4Be8FermiFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4Be8FermiFragment::operator==(const G4Be8FermiFragment &) const
{
    return false;
}

G4bool G4Be8FermiFragment::operator!=(const G4Be8FermiFragment &) const
{
    return true;
}



G4Be8FermiFragment::G4Be8FermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE)
  : G4UnstableFermiFragment(anA,aZ,Pol,ExE)
{
  // Be8 ----> alpha + alpha 
  G4double alpha_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4);

  Masses.push_back(alpha_mass);
  Masses.push_back(alpha_mass);
  
  AtomNum.push_back(4);
  AtomNum.push_back(4);
  
  Charges.push_back(2);
  Charges.push_back(2);
 
}
