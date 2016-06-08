// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Colour.cc,v 1.2.8.1 1999/12/07 20:48:53 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// John Allison 20th October 1996

#include "G4Colour.hh"
#include "G4ios.hh"

ostream& operator << (ostream& os, const G4Colour& c) {
  return os << '(' << c.red << ',' << c.green << ',' << c.blue
	    << ',' << c.alpha << ')';
}

G4bool G4Colour::operator != (const G4Colour& c) const {
  if (
      (red   != c.red)   ||
      (green != c.green) ||
      (blue  != c.blue)  ||
      (alpha != c.alpha)
      )
    return true;
  return false;
}
