// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara


#ifndef G4PreCompoundParameters_h
#define G4PreCompoundParameters_h 1

#include "globals.hh"

class G4PreCompoundParameters
{
private:
  static G4PreCompoundParameters  thePreCompoundParameters;

  // Level density parameter
  const G4double theLevelDensity;


  // default constructor
  G4PreCompoundParameters() : theLevelDensity(0.125) {}
  //  G4PreCompoundParameters(G4int Dummy) {G4int i = Dummy;}

public:

  ~G4PreCompoundParameters() {};

  static G4PreCompoundParameters * GetAddress();

  G4double GetLevelDensity()
  { return theLevelDensity; }
  
};

#endif
