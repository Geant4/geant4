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



#ifndef G4EvaporationLevelDensityParameter_h
#define G4EvaporationLevelDensityParameter_h 1


#include "G4VLevelDensityParameter.hh"

class G4EvaporationLevelDensityParameter : public G4VLevelDensityParameter
{
public:
  G4EvaporationLevelDensityParameter() :  EvapLevelDensityParameter(0.125*(1./MeV)) {};
  virtual ~G4EvaporationLevelDensityParameter() {};

private:  
  G4EvaporationLevelDensityParameter(const G4EvaporationLevelDensityParameter &right);

  const G4EvaporationLevelDensityParameter & operator=(const G4EvaporationLevelDensityParameter &right);
  G4bool operator==(const G4EvaporationLevelDensityParameter &right) const;
  G4bool operator!=(const G4EvaporationLevelDensityParameter &right) const;
  
public:
  G4double LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const 
    {return EvapLevelDensityParameter;}

private:

  const G4double EvapLevelDensityParameter;

};


#endif
