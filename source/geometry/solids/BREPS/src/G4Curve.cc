// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.cc,v 1.2 1999-12-15 14:50:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4Curve.hh"

G4Curve::G4Curve():bounded(false),bBox(G4BoundingBox3D::space),
                                  sameSense(true){}

G4Curve::~G4Curve(){}

