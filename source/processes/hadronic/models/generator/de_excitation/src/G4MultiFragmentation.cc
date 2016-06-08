// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MultiFragmentation.cc,v 1.1 1999/01/07 16:11:55 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)

#include "G4MultiFragmentation.hh"

G4MultiFragmentation::G4MultiFragmentation()
{
}

G4MultiFragmentation::G4MultiFragmentation(const G4MultiFragmentation &right)
{
}


G4MultiFragmentation::~G4MultiFragmentation()
{
}


const G4MultiFragmentation & G4MultiFragmentation::operator=(const G4MultiFragmentation &right)
{
   G4Exception("G4MultiFragmentation::operator= meant to not be accessable");
  return *this;
}


int G4MultiFragmentation::operator==(const G4MultiFragmentation &right) const
{
  return (this == (G4MultiFragmentation *) &right);
}

int G4MultiFragmentation::operator!=(const G4MultiFragmentation &right) const
{
  return (this != (G4MultiFragmentation *) &right);
}

