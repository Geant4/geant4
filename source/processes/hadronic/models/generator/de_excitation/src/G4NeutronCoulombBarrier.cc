// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NeutronCoulombBarrier.cc,v 1.1 2000/06/09 11:43:36 larazb Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4NeutronCoulombBarrier.hh"

G4NeutronCoulombBarrier::G4NeutronCoulombBarrier(const G4NeutronCoulombBarrier & right)
{
  G4Exception("G4NeutronCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4NeutronCoulombBarrier & G4NeutronCoulombBarrier::operator=(const G4NeutronCoulombBarrier & right)
{
 G4Exception("G4NeutronCoulombBarrier::operator= meant to not be accessable.");
 return *this;
}

G4bool G4NeutronCoulombBarrier::operator==(const G4NeutronCoulombBarrier & right) const 
{
 return false;
}

G4bool G4NeutronCoulombBarrier::operator!=(const G4NeutronCoulombBarrier & right) const 
{
 return true;
}

