//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
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

G4Text::G4Text (const G4Text& text):
G4VMarker (text),
fText     (text.fText),
fLayout   (text.fLayout),
fXOffset  (text.fXOffset),
fYOffset  (text.fYOffset)
{}

G4Text::~G4Text () {}

G4Text& G4Text::operator= (const G4Text& rhs)
{
  if (&rhs == this) return *this;
  G4VMarker::operator=(rhs);
  fText = rhs.fText;
  fLayout = rhs.fLayout;
  fXOffset = rhs.fXOffset;
  fYOffset = rhs.fYOffset;
  return *this;
}

std::ostream& operator<< (std::ostream& os, const G4Text& text)
{
  os << "G4Text: \"" << text.GetText()
     << "\"\n  layout " << text.GetLayout()
     << ", offset (" << text.GetXOffset() << ',' << text.GetYOffset() << ")\n"
     << (const G4VMarker&)text;
  return os;
}

std::ostream& operator<< (std::ostream& os, G4Text::Layout layout)
{
  if (layout == G4Text::left) os << "left";
  if (layout == G4Text::centre) os << "centre";
  if (layout == G4Text::right) os << "right";
  return os;
}
