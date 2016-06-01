// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4VFermiBreakUp_h
#define G4VFermiBreakUp_h 1

#include "globals.hh"
#include "G4FragmentVector.hh"

class G4VFermiBreakUp 
{
public:
  G4VFermiBreakUp();
  virtual ~G4VFermiBreakUp();
  
private:
  G4VFermiBreakUp(const G4VFermiBreakUp &right);
  
  const G4VFermiBreakUp & operator=(const G4VFermiBreakUp &right);
  G4bool operator==(const G4VFermiBreakUp &right) const;
  G4bool operator!=(const G4VFermiBreakUp &right) const;
  
public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;
};


#endif


