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
// $Id: G4Scale.cc,v 1.8 2005/07/05 14:04:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// John Allison  21st July 2001

#include "G4Scale.hh"

G4Scale::G4Scale (G4double length, const G4String& annotation,
		  G4Scale::Direction direction,
		  G4bool autoPlacing,
		  G4double xmid, G4double ymid, G4double zmid):
  fLength(length),
  fAnnotation(annotation),
  fDirection(direction),
  fAutoPlacing(autoPlacing),
  fXmid(xmid), fYmid(ymid), fZmid(zmid) {}

G4Scale::~G4Scale () {}

const G4String G4Scale::GuidanceString
(
 "An annotated line in the specified direction with tick marks at the"
 "\nend.  If autoPlacing is true it is required to be centred at the"
 "\nfront, right, bottom corner of the world space, comfortably outside"
 "\nthe existing bounding box/sphere so that existing objects do not"
 "\nobscure it.  Otherwise it is required to be drawn with mid-point at"
 "\n(xmid, ymid, zmid)."
 "\n"
 "\nThe auto placing algorithm might be:"
 "\n  x = xmin + (1 + comfort) * (xmax - xmin);"
 "\n  y = ymin - comfort * (ymax - ymin);"
 "\n  z = zmin + (1 + comfort) * (zmax - zmin);"
 "\n  if direction == x then (x - length,y,z) to (x,y,z);"
 "\n  if direction == y then (x,y,z) to (x,y + length,z);"
 "\n  if direction == z then (x,y,z - length) to (x,y,z);"
);

const G4String& G4Scale::GetGuidanceString()
{
  return GuidanceString;
}
