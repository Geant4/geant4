// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QuasmonVector.hh,v 1.7 2001-11-12 15:08:59 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4QuasmonVector_h
#define G4QuasmonVector_h 1

// ----------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QuasmonVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for a Vector of Quasmons - output of CHIPS model
// ----------------------------------------------------------------

#include "G4Quasmon.hh"
#include "g4std/vector"

typedef G4std::vector<G4Quasmon *> G4QuasmonVector;
struct DeleteQuasmon{ void operator()(G4Quasmon *aN){delete aN;} };

#endif
