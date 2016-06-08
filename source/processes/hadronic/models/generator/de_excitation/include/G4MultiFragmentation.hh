// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MultiFragmentation.hh,v 1.1 1999/01/07 16:11:44 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)

#ifndef G4MultiFragmentation_h
#define G4MultiFragmentation_h 1

#include "G4FragmentVector.hh"

class G4MultiFragmentation 
{
public:
  
  G4MultiFragmentation();
  virtual ~G4MultiFragmentation();

private:
  G4MultiFragmentation(const G4MultiFragmentation &right);

  const G4MultiFragmentation & operator=(const G4MultiFragmentation &right);
  int operator==(const G4MultiFragmentation &right) const;
  int operator!=(const G4MultiFragmentation &right) const;

public:
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;
};


#endif


