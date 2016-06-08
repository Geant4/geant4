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



#ifndef G4FissionProbability_h
#define G4FissionProbability_h 1


#include "G4VEmissionProbability.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4FissionLevelDensityParameter.hh"

class G4FissionProbability : public G4VEmissionProbability
{
public:
  // Default constructor
  G4FissionProbability() {};

  ~G4FissionProbability() {};  

private:  
  // Copy constructor
  G4FissionProbability(const G4FissionProbability &right);

  const G4FissionProbability & operator=(const G4FissionProbability &right);
  G4bool operator==(const G4FissionProbability &right) const;
  G4bool operator!=(const G4FissionProbability &right) const;
  
public:
	G4double EmissionProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy);

private:
	G4double EvaporationPairingCorrection(const G4int A, const G4int Z) const;

	G4double FissionPairingCorrection(const G4int A, const G4int Z) const;


	G4EvaporationLevelDensityParameter theEvapLDP;
	G4FissionLevelDensityParameter theFissLDP;


};


#endif
