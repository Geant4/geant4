// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisPrim.cc,v 1.4 1999-05-25 09:10:26 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  August 1995

#include "G4VVisPrim.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4VVisPrim::~G4VVisPrim () {}

G4Visible& G4VVisPrim::operator = (const G4Visible& right) {
  return G4Visible::operator = (right);
}

G4VVisPrim& G4VVisPrim::operator = (const G4VVisPrim& right) {
  if (&right == this) return *this;
  G4Visible::operator = (right);
  return *this;
}

G4bool G4VVisPrim::operator == (const G4Visible& right) const{
  return G4Visible::operator == (right);
}

G4bool G4VVisPrim::operator == (const G4VVisPrim& right) const{
  return G4Visible::operator == (right);
}

ostream& operator << (ostream& os, const G4VVisPrim& prim) {
  return os << (G4Visible) prim;
}
