// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#ifndef G4PreCompoundTransitions_h
#define G4PreCompoundTransitions_h 1

// Compute transition probailities:
// TransitionProb1 => probability of transition with  \Delta N = +2
//                    number of excitons will be increased on 2
// TransitionProb2 => probability of transition with  \Delta N = -2
//                    number of excitons will be decreased on 2
// TransitionProb3 => probability of transition with  \Delta N = 0
//                    number of excitons will be the same

#include "globals.hh"
#include "G4Fragment.hh"
#include "G4PreCompoundParameters.hh"
#include "G4Proton.hh"
#include "Randomize.hh"

class G4PreCompoundTransitions
{
public:

  // Calculates transition probabilities with Delta N = +2 (Trans1) -2 (Trans2) and 0 (Trans3)
  G4PreCompoundTransitions(const G4Fragment & aFragment);

  ~G4PreCompoundTransitions() {};

private:
  G4PreCompoundTransitions() {};
  
  G4PreCompoundTransitions(const G4PreCompoundTransitions &right) {};
  
  const G4PreCompoundTransitions& operator=(const G4PreCompoundTransitions &right);

  G4bool operator==(const G4PreCompoundTransitions &right) const;
  
  G4bool operator!=(const G4PreCompoundTransitions &right) const;


public:
  G4double GetTotalProbability(void)
  { return TransitionProb1+TransitionProb2+TransitionProb3; }
  
  G4Fragment PerformTransition(const G4Fragment & aFragment);
  
private:
	
  G4double TransitionProb1;
  G4double TransitionProb2;
  G4double TransitionProb3;

};

#endif
