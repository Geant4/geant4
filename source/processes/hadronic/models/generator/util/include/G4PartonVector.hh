// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PartonVector.hh,v 1.1 1998/08/22 08:57:57 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4PartonVector_h
#define G4PartonVector_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4PartonVector ----------------
//             by Gunter Folger, May 1998.
// ------------------------------------------------------------

#include "G4Parton.hh"
#include <rw/tpordvec.h>

typedef RWTPtrOrderedVector<G4Parton> G4PartonVector;

#endif
