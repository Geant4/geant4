// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VMarker.cc,v 1.1 1999-01-07 16:09:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4VMarker.hh"

#include "G4VisAttributes.hh"

ostream& operator << (ostream& os, const G4VMarker& marker) {
  os << "G4VMarker: position: " << marker.fPosition
     << ", world size: " << marker.fWorldSize
     << ", screen size: " << marker.fScreenSize << '\n'
     << "           fill style: ";
  switch (marker.fFillStyle) {
  case G4VMarker::noFill:
    os << "no fill";
    break;
  case G4VMarker::hashed:
    os << "hashed";
    break;
  case G4VMarker::filled:
    os << "filled";
    break;
  default:
    os << "unrecognised"; break;
  }
  os << "\n           " << (G4VVisPrim) marker;
  return os;
}

G4bool operator != (const G4VMarker& m1, const G4VMarker& m2) {
  if (
      (m1.fWorldSize != m2.fWorldSize) ||
      (m1.fScreenSize != m2.fScreenSize) ||
      (m1.fFillStyle != m2.fFillStyle) ||
      !(m1.fPosition == m2.fPosition)
      )
    return true;
  return false;
}
