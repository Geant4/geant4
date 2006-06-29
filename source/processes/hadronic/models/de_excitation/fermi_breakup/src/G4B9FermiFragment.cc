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
// $Id: G4B9FermiFragment.cc,v 1.7 2006-06-29 20:12:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
