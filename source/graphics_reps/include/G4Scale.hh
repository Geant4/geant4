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
// $Id: G4Scale.hh,v 1.2 2001-07-24 21:41:41 johna Exp $
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

  virtual ~G4Scale ();

  virtual G4Visible&  operator = (const G4Visible& from);
  virtual G4VVisPrim& operator = (const G4VVisPrim& from);
  virtual G4VMarker&  operator = (const G4VMarker& from);
  virtual G4Scale&    operator = (const G4Scale& from);

  G4double        GetLength() const;
  const G4String& GetAnnotation() const;
  Direction       GetDirection() const;
  G4bool          GetAutoPlacing() const;
  G4double        GetXmid() const;
  G4double        GetYmid() const;
  G4double        GetZmid() const;

private:
  G4double fLength;
  G4String fAnnotation;
  Direction fDirection;
  G4bool fAutoPlacing;
  G4double fXmid, fYmid, fZmid;
};

#include "G4Scale.icc"

#endif
