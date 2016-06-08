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



#ifndef G4FissionProbability_h
#define G4FissionProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4VEvaporationChannel.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4FissionLevelDensityParameter.hh"

class G4FissionProbability : public G4VEmissionProbability
{
public:
  // Only available constructor
  G4FissionProbability(G4VEvaporationChannel * aChannel)
    { theChannel = aChannel; };


  ~G4FissionProbability() {};  

private:  
  // Default constructor
  G4FissionProbability() {};

  // Copy constructor
  G4FissionProbability(const G4FissionProbability &right);

  const G4FissionProbability & operator=(const G4FissionProbability &right);
  G4bool operator==(const G4FissionProbability &right) const;
  G4bool operator!=(const G4FissionProbability &right) const;
  
public:
  G4double EmissionProbability(const G4Fragment & fragment, const G4double photonExcitation);

private:
  G4VEvaporationChannel * theChannel;

  G4EvaporationLevelDensityParameter theEvapLDP;
  G4FissionLevelDensityParameter theFissLDP;


};


#endif
