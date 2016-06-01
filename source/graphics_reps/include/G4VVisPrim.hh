// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisPrim.hh,v 2.0 1998/07/02 17:30:23 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  August 1995
// Virtual base class for Visualization Primitives
// (or Visualization Representations, as they are sometimes called).

#ifndef G4VVISPRIM_HH
#define G4VVISPRIM_HH

#include "globals.hh"
#include "G4Visible.hh"

class ostream;
class G4VisAttributes;

class G4VVisPrim: public G4Visible {

  friend ostream& operator << (ostream& os, const G4VVisPrim& prim);

public:

  G4VVisPrim ();
  G4VVisPrim (const G4VVisPrim& prim);
  G4VVisPrim& operator = (const G4VVisPrim& right);
  G4VVisPrim (const G4VisAttributes* pVA);

  virtual G4bool operator == (const G4VVisPrim& prim) const;

};

#include "G4VVisPrim.icc"

#endif
