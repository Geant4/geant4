// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Text.cc,v 1.5 2001-02-03 18:29:58 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  17/11/96.

#include "G4Text.hh"

G4Text::G4Text (const G4String& text):
fText   (text),
fLayout (left),
fXOffset(0.) , fYOffset(0.)
{}

G4Text::G4Text (const G4String& text, const G4Point3D& pos):
G4VMarker (pos),
fText     (text),
fLayout   (left),
fXOffset(0.) , fYOffset(0.)
{}

G4Text::G4Text (const G4VMarker& marker):
G4VMarker (marker),
fText     ("")    ,
fLayout   (left)  ,
fXOffset(0.) , fYOffset(0.)
{}

G4Text::~G4Text () {}

G4Visible & G4Text::operator = (const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4Text::operator = (const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4VMarker & G4Text::operator = (const G4VMarker &from) {
  return G4VMarker::operator = (from);
}

G4Text & G4Text::operator = (const G4Text &from) {
  if (&from == this) return *this;
  G4VMarker::operator = (from);
  fText    = from.fText;
  fLayout  = from.fLayout;
  fXOffset = from.fXOffset;
  fYOffset = from.fYOffset;
  return *this;
}
