// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisPrim.hh,v 1.2 1999-05-12 16:10:56 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  August 1995
// Virtual base class for Visualization Primitives
// (or Visualization Representations, as they are sometimes called).

#ifndef G4VVISPRIM_HH
#define G4VVISPRIM_HH

#include "G4Visible.hh"

class ostream;
class G4VisAttributes;

class G4VVisPrim: public G4Visible {

  friend ostream& operator << (ostream& os, const G4VVisPrim& prim);

public:

  G4VVisPrim ();
  G4VVisPrim (const G4VVisPrim& prim);
  G4VVisPrim (const G4VisAttributes* pVA);

  virtual ~G4VVisPrim ();

  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
  virtual G4bool operator == (const G4VVisPrim& right) const;

};

#include "G4VVisPrim.icc"

#endif
