// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyline.cc,v 1.2 1999-05-19 08:33:50 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

#include "G4Polyline.hh"

G4Polyline::G4Polyline (const G4VVisPrim& prim):
  G4VVisPrim (prim)
{}

G4Visible & G4Polyline::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4Polyline::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

ostream& operator << (ostream& os, const G4Polyline& line) {
  os << "G4Polyline: ";
  os << '\n' << (G4VVisPrim) line;
  os << '\n' << (G4Point3DList) line;
  return os;
}
