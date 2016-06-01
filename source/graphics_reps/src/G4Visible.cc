// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Visible.cc,v 2.1 1998/07/13 16:56:24 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
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
