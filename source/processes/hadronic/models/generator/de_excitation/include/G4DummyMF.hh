// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DummyMF.hh,v 1.1 1999/01/07 16:11:39 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)

#ifndef G4DummyMF_h
#define G4DummyMF_h 1

#include "G4MultiFragmentation.hh"

class G4DummyMF : public G4MultiFragmentation
{
public:
  G4DummyMF();
  ~G4DummyMF();
  
private:
  G4DummyMF(const G4DummyMF &right);
  const G4DummyMF & operator=(const G4DummyMF &right);
  int operator==(const G4DummyMF &right) const;
  int operator!=(const G4DummyMF &right) const;

public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
  
};

#endif
