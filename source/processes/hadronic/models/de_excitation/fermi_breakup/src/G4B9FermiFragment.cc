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
// $Id: G4B9FermiFragment.cc,v 1.5 2003/12/02 12:22:03 larazb Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4B9FermiFragment.hh"
#include "G4IonTable.hh"
#include "G4HadronicException.hh"

G4B9FermiFragment::G4B9FermiFragment()
{
}

G4B9FermiFragment::G4B9FermiFragment(const G4B9FermiFragment &) : G4UnstableFermiFragment()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4B9FermiFragment::copy_constructor meant to not be accessable");
}


G4B9FermiFragment::~G4B9FermiFragment()
{
}


const G4B9FermiFragment & G4B9FermiFragment::operator=(const G4B9FermiFragment &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4B9FermiFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4B9FermiFragment::operator==(const G4B9FermiFragment &) const
{
    return false;
}

G4bool G4B9FermiFragment::operator!=(const G4B9FermiFragment &) const
{
    return true;
}



G4B9FermiFragment::G4B9FermiFragment(const G4int anA, const G4int aZ, const G4int Pol, const G4double ExE)
  : G4UnstableFermiFragment(anA,aZ,Pol,ExE)
{
  // B9 ----> alpha + alpha + proton

  G4double alpha_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(2,4);
  G4double proton_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(1,1);
  
  Masses.push_back(alpha_mass);
  Masses.push_back(alpha_mass);
  Masses.push_back(proton_mass);
  
  AtomNum.push_back(4);
  AtomNum.push_back(4);
  AtomNum.push_back(1);
  
  Charges.push_back(2);
  Charges.push_back(2);
  Charges.push_back(1);

}
