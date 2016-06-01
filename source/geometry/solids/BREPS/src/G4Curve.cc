// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.cc,v 2.3 1998/11/06 15:31:26 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4Curve.hh"

G4Curve::G4Curve():bounded(false),bBox(G4BoundingBox3D::space),
                                  sameSense(true){}

G4Curve::~G4Curve(){}

