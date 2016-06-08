// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4VFermiFragment.hh"


G4VFermiFragment::G4VFermiFragment(const G4VFermiFragment &right)
{
  G4Exception("G4VFermiFragment::copy_constructor meant to not be accessable");
}


const G4VFermiFragment & G4VFermiFragment::operator=(const G4VFermiFragment &right)
{
  G4Exception("G4VFermiFragment::operator= meant to not be accessable");
  return *this;
}


G4bool G4VFermiFragment::operator==(const G4VFermiFragment &right) const
{
  return false;
}

G4bool G4VFermiFragment::operator!=(const G4VFermiFragment &right) const
{
  return true;
}



