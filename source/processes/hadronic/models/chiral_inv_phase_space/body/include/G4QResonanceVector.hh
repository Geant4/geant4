// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QResonanceVector.hh,v 1.2 1999-12-15 14:52:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4QResonanceVector_h
#define G4QResonanceVector_h 1

// -------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QResonanceVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for a Vector of Resonances used by CHIPS model
// --------------------------------------------------------------

#include "G4QResonance.hh"
#include <rw/tpordvec.h>

typedef RWTPtrOrderedVector<G4QResonance> G4QResonanceVector;

#endif
