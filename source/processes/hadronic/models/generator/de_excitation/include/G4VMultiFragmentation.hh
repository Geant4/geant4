// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara 

#ifndef G4VMultiFragmentation_h
#define G4VMultiFragmentation_h 1

#include "globals.hh"
#include "G4Fragment.hh"

class G4VMultiFragmentation 
{
public:
  G4VMultiFragmentation();
  virtual ~G4VMultiFragmentation();
  
private:
  G4VMultiFragmentation(const G4VMultiFragmentation &right);
  
  const G4VMultiFragmentation & operator=(const G4VMultiFragmentation &right);
  G4bool operator==(const G4VMultiFragmentation &right) const;
  G4bool operator!=(const G4VMultiFragmentation &right) const;
  
public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;
};


#endif


