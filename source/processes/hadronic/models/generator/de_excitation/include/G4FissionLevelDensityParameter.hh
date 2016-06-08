// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) 
//



#ifndef G4FissionLevelDensityParameter_h
#define G4FissionLevelDensityParameter_h 1


#include "G4VLevelDensityParameter.hh"
#include "G4EvaporationLevelDensityParameter.hh"


class G4FissionLevelDensityParameter : public G4VLevelDensityParameter
{
public:
  G4FissionLevelDensityParameter() {};
  virtual ~G4FissionLevelDensityParameter() {};

private:  
  G4FissionLevelDensityParameter(const G4FissionLevelDensityParameter &right);

  const G4FissionLevelDensityParameter & operator=(const G4FissionLevelDensityParameter &right);
  G4bool operator==(const G4FissionLevelDensityParameter &right) const;
  G4bool operator!=(const G4FissionLevelDensityParameter &right) const;
  
public:
  G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const;

  
private:
  
  G4EvaporationLevelDensityParameter theEvaporationLevelDensityParameter;

};


#endif
