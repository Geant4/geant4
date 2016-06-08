// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlacedPolyhedron.cc,v 1.1.6.1 1999/12/07 20:48:55 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $

#include "G4PlacedPolyhedron.hh"

G4PlacedPolyhedron::G4PlacedPolyhedron () {}

G4PlacedPolyhedron::G4PlacedPolyhedron
(const G4Polyhedron& polyhedron, const G4Transform3D& transform):
  fPolyhedron (polyhedron), fTransform (transform) {}
