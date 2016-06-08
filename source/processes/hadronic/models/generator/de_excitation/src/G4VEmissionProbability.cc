// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4VEmissionProbability.hh"


G4VEmissionProbability::G4VEmissionProbability(const G4VEmissionProbability &right)
{
 G4Exception("G4VEmissionProbability::copy_constructor meant to not be accessable");
}




const G4VEmissionProbability & G4VEmissionProbability::operator=(const G4VEmissionProbability &right)
{
  G4Exception("G4VEmissionProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4VEmissionProbability::operator==(const G4VEmissionProbability &right) const
{
  return false;
}

G4bool G4VEmissionProbability::operator!=(const G4VEmissionProbability &right) const
{
  return true;
}





