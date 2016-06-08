// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.cc,v 1.1.10.1 1999/12/07 20:48:23 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#include "G4Curve.hh"

G4Curve::G4Curve():bounded(false),bBox(G4BoundingBox3D::space),
                                  sameSense(true){}

G4Curve::~G4Curve(){}

