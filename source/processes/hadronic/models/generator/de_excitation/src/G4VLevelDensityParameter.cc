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


#include "G4VLevelDensityParameter.hh"


G4VLevelDensityParameter::
G4VLevelDensityParameter(const G4VLevelDensityParameter &right)
{
 G4Exception("G4VLevelDensityParameter::copy_constructor meant to not be accessable");
}




const G4VLevelDensityParameter & G4VLevelDensityParameter::
operator=(const G4VLevelDensityParameter &right)
{
  G4Exception("G4VLevelDensityParameter::operator= meant to not be accessable");
  return *this;
}


G4bool G4VLevelDensityParameter::
operator==(const G4VLevelDensityParameter &right) const
{
  return false;
}

G4bool G4VLevelDensityParameter::
operator!=(const G4VLevelDensityParameter &right) const
{
  return true;
}





