// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedStringDecay.cc,v 1.1.10.1 1999/12/07 20:51:56 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// G4ExcitedStringDecay

#include "G4ExcitedStringDecay.hh"

G4ExcitedStringDecay::G4ExcitedStringDecay()
{
	theStringDecay=NULL;
}

G4ExcitedStringDecay::G4ExcitedStringDecay(G4VLongitudinalStringDecay * aStringDecay)
: theStringDecay(aStringDecay)
{}

G4ExcitedStringDecay::G4ExcitedStringDecay(const G4ExcitedStringDecay &right)
{
}

G4ExcitedStringDecay::~G4ExcitedStringDecay()
{
}

const G4ExcitedStringDecay & G4ExcitedStringDecay::operator=(const G4ExcitedStringDecay &right)
{
  G4Exception("G4ExcitedStringDecay::operator= meant to not be accessable");
  return *this;
}

int G4ExcitedStringDecay::operator==(const G4ExcitedStringDecay &right) const
{
  return 0;
}

int G4ExcitedStringDecay::operator!=(const G4ExcitedStringDecay &right) const
{
  return 1;
}

