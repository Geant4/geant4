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
// $Id: G4Scale.cc 83392 2014-08-21 14:36:35Z gcosmo $
//
// 
// John Allison  21st July 2001

#include "G4Scale.hh"

G4Scale::G4Scale (G4double length, const G4String& annotation,
                  G4Scale::Direction direction,
                  G4bool autoPlacing,
                  G4double xmid, G4double ymid, G4double zmid,
                  G4double annotationSize):
  fLength(length),
  fAnnotation(annotation),
  fAnnotationSize(annotationSize),
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
