// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) written from G4Evaporation.hh (May 1998)
//



#ifndef G4VEvaporation_h
#define G4VEvaporation_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VEvaporation 
{
public:
  G4VEvaporation() {};
  virtual ~G4VEvaporation() {}; // *

private:  
  G4VEvaporation(const G4VEvaporation &right);

  const G4VEvaporation & operator=(const G4VEvaporation &right);
  G4bool operator==(const G4VEvaporation &right) const;
  G4bool operator!=(const G4VEvaporation &right) const;
  
public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;


};


#endif
