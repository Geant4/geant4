// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyline.hh,v 1.3 1999-05-25 09:10:09 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  virtual G4Visible&  operator = (const G4Visible& right);
  virtual G4VVisPrim& operator = (const G4VVisPrim& right);
  virtual G4Polyline& operator = (const G4Polyline& right);
};

inline G4Polyline::G4Polyline () {}

#endif
