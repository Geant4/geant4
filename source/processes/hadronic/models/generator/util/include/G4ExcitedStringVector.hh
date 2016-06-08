// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedStringVector.hh,v 1.3 1999/12/15 14:52:50 gunter Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
#ifndef G4ExcitedStringVector_h
#define G4ExcitedStringPartonVector_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4ExcitedStringVector ----------------
//             by Gunter Folger, June 1998.
// ------------------------------------------------------------

#include "G4ExcitedString.hh"
#include "g4rw/tpordvec.h"

typedef G4RWTPtrOrderedVector<G4ExcitedString> G4ExcitedStringVector;

#endif
