// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QCandidateVector.hh,v 1.1 1999-11-17 11:04:14 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4QCandidateVector_h
#define G4QCandidateVector_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Quasmon initiated Candidates used by CHIPS model
// -----------------------------------------------------------------

#include "G4QCandidate.hh"
#include <rw/tpordvec.h>

typedef RWTPtrOrderedVector<G4QCandidate> G4QCandidateVector;

#endif
