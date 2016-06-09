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
// $Id: G4Scale.hh,v 1.7 2006-06-29 19:05:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  21st July 2001.

#ifndef G4SCALE_HH
#define G4SCALE_HH

#include "G4VMarker.hh"
#include "globals.hh"

class G4Scale: public G4VMarker {

public: // With description

  enum Direction {x, y, z};

  G4Scale (G4double length, const G4String& annotation = "",
	   Direction direction = x,
	   G4bool autoPlacing = true,
	   G4double xmid = 0., G4double ymid = 0., G4double zmid = 0.);
  // This creates a representation of annotated line in the specified
  // direction with tick marks at the end.  If autoPlacing is true it
  // is required to be centred at the front, right, bottom corner of
  // the world space, comfortably outside the existing bounding
  // box/sphere so that existing objects do not obscure it.  Otherwise
  // it is required to be drawn with mid-point at (xmid, ymid, zmid).
  //
  // The auto placing algorithm might be:
  //   x = xmin + (1 + comfort) * (xmax - xmin)
  //   y = ymin - comfort * (ymax - ymin)
  //   z = zmin + (1 + comfort) * (zmax - zmin)
  //   if direction == x then (x - length,y,z) to (x,y,z)
  //   if direction == y then (x,y,z) to (x,y + length,z)
  //   if direction == z then (x,y,z - length) to (x,y,z)

  // Uses compiler-generated copy constructor.

  ~G4Scale ();

  G4double        GetLength() const;
  const G4String& GetAnnotation() const;
  Direction       GetDirection() const;
  G4bool          GetAutoPlacing() const;
  G4double        GetXmid() const;
  G4double        GetYmid() const;
  G4double        GetZmid() const;

  static const G4String& GetGuidanceString();

private:

  G4double fLength;
  G4String fAnnotation;
  Direction fDirection;
  G4bool fAutoPlacing;
  G4double fXmid, fYmid, fZmid;
  static const G4String GuidanceString;
};

#include "G4Scale.icc"

#endif
