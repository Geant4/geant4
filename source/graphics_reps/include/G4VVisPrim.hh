// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisPrim.hh,v 1.6 1999-12-15 16:46:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/iostream"

class G4VisAttributes;

class G4VVisPrim: public G4Visible {

  friend G4std::ostream& operator << (G4std::ostream& os, const G4VVisPrim& prim);

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
