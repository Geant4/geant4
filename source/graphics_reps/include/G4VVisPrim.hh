// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisPrim.hh,v 1.4.2.1 1999/12/07 20:48:51 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// John Allison  August 1995

// Class Description:
// Virtual base class for Visualization Primitives
// (or Visualization Representations, as they are sometimes called).
// Class Description - End:

#ifndef G4VVISPRIM_HH
#define G4VVISPRIM_HH

#include "G4Visible.hh"

class ostream;
class G4VisAttributes;

class G4VVisPrim: public G4Visible {

  friend ostream& operator << (ostream& os, const G4VVisPrim& prim);

public: // With description

  G4VVisPrim ();
  G4VVisPrim (const G4VVisPrim& prim);
  G4VVisPrim (const G4VisAttributes* pVA);

  virtual ~G4VVisPrim ();

  virtual G4Visible& operator = (const G4Visible& right);
  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
  virtual G4bool operator == (const G4Visible& right) const;
  virtual G4bool operator == (const G4VVisPrim& right) const;

};

#include "G4VVisPrim.icc"

#endif
