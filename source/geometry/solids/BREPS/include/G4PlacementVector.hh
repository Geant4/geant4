// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlacementVector.hh,v 1.3 2000-08-28 08:57:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PlacementVector
//
// Class description:
//
// A ordered collection of pointers to placement axes
// (G4Axis2Placement3D) in the 3D space.

// Authors: J.Sulkimo, P.Urban.
// ----------------------------------------------------------------------
#ifndef included_G4PlacementVector
#define included_G4PlacementVector

#include "g4rw/tpordvec.h"
#include "G4Axis2Placement3D.hh"

typedef G4RWTPtrOrderedVector<G4Axis2Placement3D> G4PlacementVector;

#endif
