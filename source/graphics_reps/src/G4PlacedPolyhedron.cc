// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlacedPolyhedron.cc,v 1.1 1999-07-21 16:43:48 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4PlacedPolyhedron.hh"

G4PlacedPolyhedron::G4PlacedPolyhedron () {}

G4PlacedPolyhedron::G4PlacedPolyhedron
(const G4Polyhedron& polyhedron, const G4Transform3D& transform):
  fPolyhedron (polyhedron), fTransform (transform) {}
