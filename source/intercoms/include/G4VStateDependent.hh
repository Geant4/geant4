// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VStateDependent.hh,v 1.1 1999-01-07 16:09:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VStateDependent ----------------
//             by Gabriele Cosmo, November 1996
// Abstract class responsible to Notify the change of state
// ------------------------------------------------------------

#ifndef G4VStateDependent_h
#define G4VStateDependent_h 1

#include "globals.hh"
#include "G4ApplicationState.hh"

class G4VStateDependent {

public:
  G4VStateDependent();
  virtual ~G4VStateDependent();
  G4int operator==(const G4VStateDependent &right) const;
  G4int operator!=(const G4VStateDependent &right) const;

  virtual G4bool Notify(G4ApplicationState requestedState) = 0;

private:
  G4VStateDependent(const G4VStateDependent &right);
  G4VStateDependent& operator=(const G4VStateDependent &right);

};

#endif
