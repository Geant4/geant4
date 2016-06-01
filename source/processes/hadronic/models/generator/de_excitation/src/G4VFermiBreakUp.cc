// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4VFermiBreakUp.hh"


G4VFermiBreakUp::G4VFermiBreakUp()
{
}

G4VFermiBreakUp::G4VFermiBreakUp(const G4VFermiBreakUp &right)
{
  G4Exception("G4VFermiBreakUp::copy_constructor meant to not be accessable");
}


G4VFermiBreakUp::~G4VFermiBreakUp()
{
}


const G4VFermiBreakUp & G4VFermiBreakUp::operator=(const G4VFermiBreakUp &right)
{
  G4Exception("G4VFermiBreakUp::operator= meant to not be accessable");
  return *this;
}


G4bool G4VFermiBreakUp::operator==(const G4VFermiBreakUp &right) const
{
  return false;
}

G4bool G4VFermiBreakUp::operator!=(const G4VFermiBreakUp &right) const
{
  return true;
}



