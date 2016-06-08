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


#ifndef G4VEvaporationChannel_h
#define G4VEvaporationChannel_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VEvaporationChannel
{
public:
  G4VEvaporationChannel() {};
  virtual ~G4VEvaporationChannel() {};

private:
  G4VEvaporationChannel(const G4VEvaporationChannel & right);

  const G4VEvaporationChannel & operator=(const G4VEvaporationChannel & right);
public:
  G4bool operator==(const G4VEvaporationChannel & right) const;
  G4bool operator!=(const G4VEvaporationChannel & right) const;

public:
  virtual void Initialize(const G4Fragment & fragment) = 0;

  virtual G4FragmentVector * BreakUp(const G4Fragment & theNucleus) = 0;

  virtual G4double GetEmissionProbability(void) const = 0;


};


#endif
