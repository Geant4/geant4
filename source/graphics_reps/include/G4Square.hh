// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Square.hh,v 1.1 1999-01-07 16:09:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  17/11/96.

#ifndef G4SQUARE_HH
#define G4SQUARE_HH

#include "G4VMarker.hh"

class G4Square: public G4VMarker {
public:
  G4Square ();
  G4Square (const G4Point3D& pos);
  G4Square (const G4VMarker& marker);
};

#include "G4Square.icc"

#endif
