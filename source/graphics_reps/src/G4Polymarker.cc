// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polymarker.cc,v 1.3 1999-05-19 08:33:51 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  November 1996

#include "G4Polymarker.hh"

G4Polymarker::G4Polymarker ():
fMarkerType (line)
{}

G4Visible & G4Polymarker::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4Polymarker::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

G4VMarker & G4Polymarker::operator = (const G4VMarker &right) {
  return G4VMarker::operator = (right);
}

ostream& operator << (ostream& os, const G4Polymarker& marker) {
  os << "G4Polymarker: type: ";
  switch (marker.fMarkerType) {
  case G4Polymarker::line:
    os << "line"; break;
  case G4Polymarker::dots:
    os << "dots"; break;
  case G4Polymarker::circles:
    os << "circles"; break;
  case G4Polymarker::squares:
    os << "squares"; break;
  default:
    os << "unrecognised"; break;
  }
  os << "\n  ";
  os << (G4VMarker) marker;
  os << (G4Point3DList) marker;
  return os;
}
