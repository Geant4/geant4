// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlacedPolyhedron.hh,v 1.3 1999-10-05 07:50:16 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4PLACEDPOLYHEDRON_HH
#define G4PLACEDPOLYHEDRON_HH

#include "G4Polyhedron.hh"
#include "G4Transform3D.hh"
#include <rw/tvordvec.h>

class G4PlacedPolyhedron {

public:

  G4PlacedPolyhedron ();
  G4PlacedPolyhedron (const G4Polyhedron&, const G4Transform3D&);

  // Uses default copy constructor, destructor and assignment.

  G4bool G4PlacedPolyhedron::operator ==
  (const G4PlacedPolyhedron& right) const {
    return this == &right;
  }

  const G4Polyhedron&  GetPolyhedron () const {return fPolyhedron;}
  const G4Transform3D& GetTransform  () const {return fTransform;}

  void SetPolyhedron (const G4Polyhedron& polyhedron) {
    fPolyhedron = polyhedron;
  }
  void SetTransform  (const G4Transform3D& transform) {
    fTransform = transform;
  }

private:

  G4Polyhedron fPolyhedron;
  G4Transform3D fTransform;

};

typedef RWTValOrderedVector<G4PlacedPolyhedron> G4PlacedPolyhedronList;

#endif /* G4PLACEDPOLYHEDRON_HH */
