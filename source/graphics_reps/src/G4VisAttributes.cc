// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisAttributes.cc,v 1.1 1999-01-07 16:09:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  23rd October 1996

#include "G4VisAttributes.hh"

G4VisAttributes::G4VisAttributes ():
fVisible           (true),
fColour            (G4Colour ()),
fLineStyle         (unbroken),
fLineWidth         (1.),
fForceDrawingStyle (false)
{}

G4VisAttributes::G4VisAttributes (G4bool visibility):
fVisible           (visibility),
fColour            (G4Colour ()),
fLineStyle         (unbroken),
fLineWidth         (1.),
fForceDrawingStyle (false)
{}

G4VisAttributes::G4VisAttributes (const G4Colour& colour):
fVisible           (true),
fColour            (colour),
fLineStyle         (unbroken),
fLineWidth         (1.),
fForceDrawingStyle (false)
{}

G4VisAttributes::G4VisAttributes (G4bool visibility,
					 const G4Colour& colour):
fVisible           (visibility),
fColour            (colour),
fLineStyle         (unbroken),
fLineWidth         (1.),
fForceDrawingStyle (false)
{}

const G4VisAttributes  G4VisAttributes::Invisible = G4VisAttributes (false);

ostream& operator << (ostream& os, const G4VisAttributes& a) {
  
  os << "G4VisAttributes: ";
  if (&a){
    if (!a.fVisible) os << " in";
    os << "visible, colour: " << a.fColour;
    os << "\n  linestyle: ";
    switch (a.fLineStyle) {
    case G4VisAttributes::unbroken:
      os << "solid"; break;
    case G4VisAttributes::dashed:
      os << "dashed"; break;
    case G4VisAttributes::dotted: os << "dotted"; break;
    default: os << "unrecognised"; break;
    }
    os << ", line width: " << a.fLineWidth;
    os << "\n  drawing style ";
    if (a.fForceDrawingStyle) {
      os << "forced to: ";
      switch (a.fForcedStyle) {
      case G4VisAttributes::wireframe:
	os << "wireframe"; break;
      case G4VisAttributes::solid:
	os << "solid"; break;
      default: os << "unrecognised"; break;
      }
    }
    else {
      os << "unforced";
    }
  } 
  else os << " The pointer is zero ";
  return os;

}

G4bool operator != (const G4VisAttributes& a1,
		    const G4VisAttributes& a2) {

  if (
      (a1.fVisible            != a2.fVisible)            ||
      (a1.fColour             != a2.fColour)             ||
      (a1.fLineStyle          != a2.fLineStyle)          ||
      (a1.fLineWidth          != a2.fLineWidth)          ||
      (a1.fForceDrawingStyle  != a2.fForceDrawingStyle)
      )
    return true;

  if (a1.fForceDrawingStyle) {
    if (a1.fForcedStyle != a2.fForcedStyle) return true;
  }

  return false;
}
