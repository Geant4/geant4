// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedStringDecay.cc,v 1.2 1998/12/09 07:53:44 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
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

