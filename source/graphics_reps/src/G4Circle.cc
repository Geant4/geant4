// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Circle.cc,v 1.2 1999-12-15 14:50:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4Circle.hh"

#include "G4VisAttributes.hh"

G4Circle::~G4Circle () {}

G4Visible & G4Circle::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4Circle::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

G4VMarker & G4Circle::operator = (const G4VMarker &right) {
  return G4VMarker::operator = (right);
}

G4Circle& G4Circle::operator = (const G4Circle& right) {
  if (&right == this) return *this;
  G4VMarker::operator = (right);
  return *this;
}
