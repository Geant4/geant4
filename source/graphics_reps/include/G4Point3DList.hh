// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DList.hh,v 1.6 1999-12-15 16:46:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

// Class Description:
// A set of 3D points.
// Class Description - End:

#ifndef G4POINT3DLIST_HH
#define G4POINT3DLIST_HH

#include "g4rw/tvordvec.h"
#include "G4Point3D.hh"
#include "g4std/iostream"

class G4Point3DList: public G4RWTValOrderedVector<G4Point3D> {

friend G4std::ostream& operator << (G4std::ostream& os, const G4Point3DList& points);

public:

  virtual ~G4Point3DList();

};

#endif
