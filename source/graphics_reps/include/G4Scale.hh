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
// $Id: G4Scale.hh,v 1.1 2001-07-22 00:49:35 johna Exp $
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
	   G4double x0 = 0, G4double y0 = 0, G4double z0 = 0);
  // This creates a representation of annotated line in the specified
  // direction with tick marks at the end.  If autoPlacing is true it
  // is required to be centred at the front, right, bottom corner of
  // the world space, say, one tenth of the way in, i.e.,
  // 
  //   for a margin, of, say 1%, so that it is just inside viewing volume:
  //   x = xmin + (1 - margin) * (xmax - xmin)
  //   y = ymin + margin * (ymax - ymin)
  //   z = zmin + (1 - margin) * (zmax - zmin)
  //   if direction == "x" then (x - length,y,z) to (x,y,z)
  //   if direction == "y" then (x,y,z) to (x,y + length,z)
  //   if direction == "z" then (x,y,z - length) to (x,y,z)
  //
  // Otherwise it is required to be drawn at (x0, y0, z0).

  G4Scale(const G4Scale&);

  virtual ~G4Scale ();

  virtual G4Visible&  operator = (const G4Visible& from);
  virtual G4VVisPrim& operator = (const G4VVisPrim& from);
  virtual G4VMarker&  operator = (const G4VMarker& from);
  virtual G4Scale&    operator = (const G4Scale& from);

  G4double        GetLength() const;
  const G4String& GetAnnotation() const;
  Direction       GetDirection() const;
  G4bool          GetAutoPlacing() const;
  G4double        GetX() const;
  G4double        GetY() const;
  G4double        GetZ() const;

private:
  G4double fLength;
  G4String fAnnotation;
  Direction fDirection;
  G4bool fAutoPlacing;
  G4double fX0, fY0, fZ0;
};

#include "G4Scale.icc"

#endif
