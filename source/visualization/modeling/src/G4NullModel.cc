// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NullModel.cc,v 1.2 1999-01-11 00:48:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  4th April 1998.
// Null model - simply a holder for modeling parameter.
// DO NOT INVOKE DescribeYourself.

#include "G4NullModel.hh"

G4NullModel::G4NullModel (const G4ModelingParameters* pMP):
  G4VModel (G4Transform3D::Identity, pMP) {}

G4NullModel::~G4NullModel () {}

void G4NullModel::DescribeYourselfTo (G4VGraphicsScene& scene) {
  G4Exception ("G4NullModel::DescribeYourselfTo called.");
}

G4bool G4NullModel::Validate () {
  return true;
}
