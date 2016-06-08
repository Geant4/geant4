// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4VLevelDensityParameter_h
#define G4VLevelDensityParameter_h 1


#include "globals.hh"

class G4VLevelDensityParameter 
{
public:
  G4VLevelDensityParameter() {};
  virtual ~G4VLevelDensityParameter() {};

private:  
  G4VLevelDensityParameter(const G4VLevelDensityParameter &right);

  const G4VLevelDensityParameter & operator=(const G4VLevelDensityParameter &right);
  G4bool operator==(const G4VLevelDensityParameter &right) const;
  G4bool operator!=(const G4VLevelDensityParameter &right) const;
  
public:
  virtual G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const = 0;

};


#endif
