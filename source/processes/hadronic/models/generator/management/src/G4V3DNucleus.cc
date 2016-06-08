// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4V3DNucleus.cc,v 1.1.10.1 1999/12/07 20:51:45 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#include "G4V3DNucleus.hh"

G4V3DNucleus::G4V3DNucleus()
{
}

G4V3DNucleus::G4V3DNucleus(const G4V3DNucleus &right)
{
}


G4V3DNucleus::~G4V3DNucleus()
{
}


const G4V3DNucleus & G4V3DNucleus::operator=(const G4V3DNucleus &right)
{
  G4Exception("G4V3DNucleus::operator= meant to not be accessable"); // needs to be looked at @@
  return *this;
}


int G4V3DNucleus::operator==(const G4V3DNucleus &right) const
{
  return 0;
}

int G4V3DNucleus::operator!=(const G4V3DNucleus &right) const
{
  return 1;
}


