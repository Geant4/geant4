// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RotationMatrix.hh,v 1.2 1999-11-16 17:40:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
