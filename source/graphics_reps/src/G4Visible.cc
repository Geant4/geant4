// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Visible.cc,v 1.1 1999-01-07 16:09:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  30th October 1996
// Base class for all things visible, i.e., which have Vis Attributes.

#include "G4Visible.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

ostream& operator << (ostream& os, const G4Visible& v) {
  if (v.fpVisAttributes) return os << *(v.fpVisAttributes);
  else return os << "No Visualization Attributes";
}
