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
// $Id: G4Text.cc,v 1.8 2006-03-13 12:46:23 allison Exp $
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

std::ostream& operator<< (std::ostream& os, const G4Text& text)
{
  os << "G4Text: \"" << text.GetText()
     << "\"\n  layout " << text.GetLayout()
     << ", offset (" << text.GetXOffset() << ',' << text.GetYOffset() << ")\n"
     << (const G4VMarker&)text;
  return os;
}
