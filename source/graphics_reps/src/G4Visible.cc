// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Visible.cc,v 1.2 1999-05-12 16:11:08 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  30th October 1996
// Base class for all things visible, i.e., which have Vis Attributes.

#include "G4Visible.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Visible::~G4Visible () {}

G4Visible& G4Visible::operator = (const G4Visible& right) {
  if (&right == this) return *this;
  fpVisAttributes = right.fpVisAttributes;
  return *this;
}

G4bool G4Visible::operator == (const G4Visible& right) const{
  return fpVisAttributes == right.fpVisAttributes;
}

ostream& operator << (ostream& os, const G4Visible& v) {
  if (v.fpVisAttributes) return os << *(v.fpVisAttributes);
  else return os << "No Visualization Attributes";
}
