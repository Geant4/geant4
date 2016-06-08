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



#ifndef G4VEmissionProbability_h
#define G4VEmissionProbability_h 1


#include "globals.hh"
#include "G4Fragment.hh"

class G4VEmissionProbability 
{
public:
  G4VEmissionProbability() {};
  virtual ~G4VEmissionProbability() {};

private:  
  G4VEmissionProbability(const G4VEmissionProbability &right);

  const G4VEmissionProbability & operator=(const G4VEmissionProbability &right);
  G4bool operator==(const G4VEmissionProbability &right) const;
  G4bool operator!=(const G4VEmissionProbability &right) const;
  
public:
  virtual G4double EmissionProbability(const G4Fragment & fragment, const G4double anEnergy) = 0;

};


#endif
