// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Visible.hh,v 1.1 1999-01-07 16:09:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  30th October 1996
// Base class for all things visible, i.e., which have Vis Attributes.

// Note: a null pointer implies no attributes.  Under those circumstances
// the visualization system is free to choose some.

#ifndef G4VISIBLE_HH
#define G4VISIBLE_HH

class G4VisAttributes;
class ostream;

class G4Visible {

  friend ostream& operator << (ostream& os, const G4Visible& v);

public:

  G4Visible ();
  G4Visible (const G4Visible& visible);
  G4Visible& operator = (const G4Visible& right);
  G4Visible (const G4VisAttributes* pVA);

  const G4VisAttributes* GetVisAttributes () const;

  void SetVisAttributes (const G4VisAttributes* pVA);
  void SetVisAttributes (const G4VisAttributes& VA);

protected:

  const G4VisAttributes* fpVisAttributes;

};

#include "G4Visible.icc"

#endif
