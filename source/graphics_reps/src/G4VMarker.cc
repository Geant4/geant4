// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VMarker.cc,v 1.3 1999-05-19 08:33:52 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4VMarker.hh"

#include "G4VisAttributes.hh"

G4VMarker::~G4VMarker () {}

G4Visible & G4VMarker::operator = (const G4Visible &right) {
  return G4Visible::operator = (right);
}

G4VVisPrim & G4VMarker::operator = (const G4VVisPrim &right) {
  return G4VVisPrim::operator = (right);
}

G4VMarker& G4VMarker::operator = (const G4VMarker& right) {
  if (&right == this) return *this;
  G4VVisPrim::operator = (right);
  fPosition   = right.fPosition;
  fWorldSize  = right.fWorldSize;  
  fScreenSize = right.fScreenSize;
  fFillStyle  = right.fFillStyle;
  fInfo       = right.fInfo;
  return *this;
}

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
      !((G4VVisPrim) m1 == (G4VVisPrim) m2)  ||
      (m1.fWorldSize != m2.fWorldSize)       ||
      (m1.fScreenSize != m2.fScreenSize)     ||
      (m1.fFillStyle != m2.fFillStyle)       ||
      !(m1.fPosition == m2.fPosition)
      )
    return true;
  return false;
}

const G4String& G4VMarker::GetInfo() const { return fInfo ;}

void G4VMarker::SetInfo( const G4String& info ){ fInfo = info ;}

