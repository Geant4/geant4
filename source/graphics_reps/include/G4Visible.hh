//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Visible.hh,v 1.7 2001-07-11 10:01:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  30th October 1996

// Class Description:
// Base class for all things visible, i.e., which have Vis Attributes.
//
// Note: a null pointer implies no attributes.  Under those circumstances
// the visualization system is free to choose some.
// Class Description - End:


#ifndef G4VISIBLE_HH
#define G4VISIBLE_HH

#include "globals.hh"
#include "g4std/iostream"

class G4VisAttributes;

class G4Visible {

  friend G4std::ostream& operator << (G4std::ostream& os, const G4Visible& v);

public: // With description

  G4Visible ();
  G4Visible (const G4Visible& visible);
  G4Visible (const G4VisAttributes* pVA);

  virtual ~G4Visible ();

  virtual G4Visible& operator = (const G4Visible& right);
  virtual G4bool operator == (const G4Visible& right) const;

  const G4VisAttributes* GetVisAttributes () const;

  void SetVisAttributes (const G4VisAttributes* pVA);
  void SetVisAttributes (const G4VisAttributes& VA);
  // The G4VisAttributes object is not stored in a G4Visible; only a
  // reference, a const pointer, is kept.  Therefore the
  // G4VisAttributes object to which it refers must have a life long
  // enough to satisfy all uses of the G4Visible object.  E.g., if the
  // G4Visible object is created on the heap (using `new') then the
  // associated G4VisAttributes object would normally also be created
  // on the heap and managed in the same way.

protected:

  const G4VisAttributes* fpVisAttributes;

};

#include "G4Visible.icc"

#endif
