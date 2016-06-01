// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
//    Constant level density parameter (for photon evaporation)
//
// by C. Dallapiccola (Nov 1998)
//


#include "G4ConstantLevelDensityParameter.hh"

G4ConstantLevelDensityParameter::
G4ConstantLevelDensityParameter(const G4ConstantLevelDensityParameter& right) :
  EvapLevelDensityParameter(0.125*(1./MeV))
{
  G4Exception("G4ConstantLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4ConstantLevelDensityParameter & G4ConstantLevelDensityParameter::
operator=(const G4ConstantLevelDensityParameter &right)
{
  G4Exception("G4ConstantLevelDensityParameter::operator= meant to not be accessable");
  return *this;
}


G4bool G4ConstantLevelDensityParameter::operator==(const G4ConstantLevelDensityParameter &right) const
{
  return false;
}

G4bool G4ConstantLevelDensityParameter::operator!=(const G4ConstantLevelDensityParameter &right) const
{
  return true;
}





