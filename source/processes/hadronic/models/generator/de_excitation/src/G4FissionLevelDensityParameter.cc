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


#include "G4FissionLevelDensityParameter.hh"


G4FissionLevelDensityParameter::
G4FissionLevelDensityParameter(const G4FissionLevelDensityParameter &right)
{
 G4Exception("G4FissionLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4FissionLevelDensityParameter & G4FissionLevelDensityParameter::
operator=(const G4FissionLevelDensityParameter &right)
{
  G4Exception("G4FissionLevelDensityParameter::operator= meant to not be accessable");
  return *this;
}


G4bool G4FissionLevelDensityParameter::
operator==(const G4FissionLevelDensityParameter &right) const
{
  return false;
}

G4bool G4FissionLevelDensityParameter::
operator!=(const G4FissionLevelDensityParameter &right) const
{
  return true;
}


G4double G4FissionLevelDensityParameter::
LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const 
{
  G4double EvapLDP = theEvaporationLevelDensityParameter.LevelDensityParameter(A,Z,U);

  if (Z >= 89) return 1.04*EvapLDP;
  else if (Z >= 85) return (1.04*(1./MeV) + 0.01*(89-Z))*EvapLDP;
  else return 1.08*EvapLDP;
}
