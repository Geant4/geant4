// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Circle.hh,v 1.1 1999-01-07 16:09:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  17/11/96.

#ifndef G4CIRCLE_HH
#define G4CIRCLE_HH

#include "G4VMarker.hh"

class G4Circle: public G4VMarker {
public:
  G4Circle ();
  G4Circle (const G4Point3D& pos);
  G4Circle (const G4VMarker& marker);
};

#include "G4Circle.icc"

#endif
