// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Colour.cc,v 1.1 1999-01-07 16:09:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison 20th October 1996

#include "G4Colour.hh"
#include "G4ios.hh"

ostream& operator << (ostream& os, const G4Colour& c) {
  return os << '(' << c.red << ',' << c.green << ',' << c.blue
	    << ',' << c.alpha << ')';
}

G4bool operator != (const G4Colour& c1, const G4Colour& c2) {
  if (
      (c1.red   != c2.red)   ||
      (c1.green != c2.green) ||
      (c1.blue  != c2.blue)  ||
      (c1.alpha != c2.alpha)
      )
    return true;
  return false;
}
