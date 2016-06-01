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


#include "G4EvaporationLevelDensityParameter.hh"


//const G4double G4EvaporationLevelDensityParameter::EvapLevelDensityParameter = 0.125*(1./MeV);

G4EvaporationLevelDensityParameter::
G4EvaporationLevelDensityParameter(const G4EvaporationLevelDensityParameter &right) :
  EvapLevelDensityParameter(0.125*(1./MeV))
{
  G4Exception("G4EvaporationLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4EvaporationLevelDensityParameter & G4EvaporationLevelDensityParameter::
operator=(const G4EvaporationLevelDensityParameter &right)
{
  G4Exception("G4EvaporationLevelDensityParameter::operator= meant to not be accessable");
  return *this;
}


G4bool G4EvaporationLevelDensityParameter::operator==(const G4EvaporationLevelDensityParameter &right) const
{
  return false;
}

G4bool G4EvaporationLevelDensityParameter::operator!=(const G4EvaporationLevelDensityParameter &right) const
{
  return true;
}





