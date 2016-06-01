// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyline.hh,v 2.0 1998/07/02 17:30:06 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  July 1995

#ifndef G4POLYLINE_HH
#define G4POLYLINE_HH

#include "G4VVisPrim.hh"
#include "G4Point3DList.hh"

class G4Polyline: public G4VVisPrim, public G4Point3DList {
  friend ostream& operator << (ostream& os, const G4Polyline& line);
public:
  G4Polyline ();
  G4Polyline (const G4VVisPrim& prim);
};

inline G4Polyline::G4Polyline () {}

#endif
