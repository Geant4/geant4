// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Point3DVector.hh,v 1.4 2000-08-28 08:57:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Point3DVector
//
// Class description:
//
// A value collection of points in 3D space (G4Point3D).

// Authors: J.Sulkimo, P.Urban.
// ----------------------------------------------------------------------
#ifndef included_G4Point3DVector
#define included_G4Point3DVector

#include "g4rw/tvvector.h"
#include "G4Point3D.hh"

typedef G4RWTValVector<G4Point3D> G4Point3DVector;

#endif
