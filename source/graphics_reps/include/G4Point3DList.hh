// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DList.hh,v 1.1 1999-01-07 16:09:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

#ifndef G4POINT3DLIST_HH
#define G4POINT3DLIST_HH

#include <rw/tvordvec.h>
#include "G4Point3D.hh"

class ostream;

class G4Point3DList: public RWTValOrderedVector<G4Point3D> {

friend ostream& operator << (ostream& os, const G4Point3DList& points);

public:

  virtual ~G4Point3DList();
  //  Destructor.

};

#include "G4Point3DList.icc"

#endif
