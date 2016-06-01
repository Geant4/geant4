// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RotationMatrix.hh,v 2.0 1998/07/02 17:33:01 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// ----------------------------------------------------------------------
//
// G4RotationMatrix class, typedef to CLHEP HepRotation
//
// ----------------------------------------------------------------------

#ifndef G4ROTATIONMATRIX_HH
#define G4ROTATIONMATRIX_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <CLHEP/Vector/Rotation.h>

typedef HepRotation G4RotationMatrix;

#endif
