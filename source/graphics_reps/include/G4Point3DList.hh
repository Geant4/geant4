// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DList.hh,v 1.7 2001-05-03 11:06:09 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  July 1995

// Class Description:
// A set of 3D points.
// Class Description - End:

#ifndef G4POINT3DLIST_HH
#define G4POINT3DLIST_HH

#include "g4std/vector"
#include "G4Point3D.hh"
#include "g4std/iostream"

class G4Point3DList: public G4std::vector<G4Point3D> {

friend G4std::ostream& operator << (G4std::ostream& os, const G4Point3DList& points);

public:

  virtual ~G4Point3DList();

};

#endif
