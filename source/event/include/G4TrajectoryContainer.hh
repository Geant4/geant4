// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoryContainer.hh,v 1.2 1999-04-15 08:41:51 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//	G4TrajectoryContainer
//

#ifndef G4TrajectoryContainer_h
#define G4TrajectoryContainer_h 1

// TrackManagement
#include "G4VTrajectory.hh"
// RWTPtrOrderedVector
#include <rw/tpordvec.h>

typedef RWTPtrOrderedVector<G4VTrajectory> G4TrajectoryContainer;

#endif


