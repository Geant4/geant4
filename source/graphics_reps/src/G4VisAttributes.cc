//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisAttributes.cc,v 1.8 2002-11-20 14:18:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  23rd October 1996

#include "G4VisAttributes.hh"

G4VisAttributes::G4VisAttributes ():
fVisible            (true),
fDaughtersInvisible (false),
fColour             (G4Colour ()),
fLineStyle          (unbroken),
fLineWidth          (1.),
fForceDrawingStyle  (false),
fAttValues          (0),
fAttDefs            (0)
{}

G4VisAttributes::G4VisAttributes (G4bool visibility):
fVisible            (visibility),
fDaughtersInvisible (false),
fColour             (G4Colour ()),
fLineStyle          (unbroken),
fLineWidth          (1.),
fForceDrawingStyle  (false),
fAttValues          (0),
fAttDefs            (0)
{}

G4VisAttributes::G4VisAttributes (const G4Colour& colour):
fVisible            (true),
fDaughtersInvisible (false),
fColour             (colour),
fLineStyle          (unbroken),
fLineWidth          (1.),
fForceDrawingStyle  (false),
fAttValues          (0),
fAttDefs            (0)
{}

G4VisAttributes::G4VisAttributes (G4bool visibility,
				  const G4Colour& colour):
fVisible            (visibility),
fDaughtersInvisible (false),
fColour             (colour),
fLineStyle          (unbroken),
fLineWidth          (1.),
fForceDrawingStyle  (false),
fAttValues          (0),
fAttDefs            (0)
{}

const G4VisAttributes  G4VisAttributes::Invisible = G4VisAttributes (false);

const G4VisAttributes& G4VisAttributes::GetInvisible() {
  return Invisible;
}

G4std::ostream& operator << (G4std::ostream& os, const G4VisAttributes& a) {
  
  os << "G4VisAttributes: ";
  if (&a){
    if (!a.fVisible) os << "in";
    os << "visible, daughters ";
    if (a.fDaughtersInvisible) os << "in";
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
    os << "\n  vector<G4AttValue> pointer is ";
    if (a.fAttValues) {
      os << "non-";
    }
    os << "zero";      
    os << "\n  vector<G4AttDef> pointer is ";
    if (a.fAttDefs) {
      os << "non-";
    }
    os << "zero";      
  } 
  else os << " zero G4VisAttributes pointer";
  return os;
}

G4bool G4VisAttributes::operator != (const G4VisAttributes& a) const {

  if (
      (fVisible            != a.fVisible)            ||
      (fDaughtersInvisible != a.fDaughtersInvisible) ||
      (fColour             != a.fColour)             ||
      (fLineStyle          != a.fLineStyle)          ||
      (fLineWidth          != a.fLineWidth)          ||
      (fForceDrawingStyle  != a.fForceDrawingStyle)  ||
      (fAttValues          != a.fAttValues)          ||
      (fAttDefs            != a.fAttDefs)
      )
    return true;

  if (fForceDrawingStyle) {
    if (fForcedStyle != a.fForcedStyle) return true;
  }

  return false;
}

G4bool G4VisAttributes::operator == (const G4VisAttributes& a) const {
  return !(G4VisAttributes::operator != (a));
}
