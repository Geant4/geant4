// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoryContainer.hh,v 1.3 1999-11-05 04:16:18 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//	G4TrajectoryContainer
//

#ifndef G4TrajectoryContainer_h
#define G4TrajectoryContainer_h 1

#include "G4VTrajectory.hh"
// RWTPtrOrderedVector
#include <rw/tpordvec.h>

// class description:
//
//  This is a container of G4VTrajectory objects and the object of this
// container will be associated to G4Event object.

typedef RWTPtrOrderedVector<G4VTrajectory> G4TrajectoryContainer;

#endif


