// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Colour.cc,v 1.3 1999-12-15 14:50:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison 20th October 1996

#include "G4Colour.hh"
#include "G4ios.hh"

G4std::ostream& operator << (G4std::ostream& os, const G4Colour& c) {
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
