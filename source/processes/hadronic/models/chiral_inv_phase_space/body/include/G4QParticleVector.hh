// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QParticleVector.hh,v 1.6 2001-10-04 20:00:21 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4QParticleVector_h
#define G4QParticleVector_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Decay Channel Vector in CHIPS model
// ---------------------------------------------------------------

#include "G4QParticle.hh"
#include "g4std/vector"

typedef G4std::vector<G4QParticle *> G4QParticleVector;
struct DeleteQParticle{ void operator()(G4QParticle *aN){delete aN;} };

#endif


