// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleVector.hh,v 1.1 1999-01-07 16:13:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
// HPW decoupling theo models from RW (Mon Mar 16 1998)
// ------------------------------------------------------------

#ifndef G4ParticleVector_h
#define G4ParticleVector_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4DynamicParticle.hh"
#include <rw/tpordvec.h>

// #ifdef STL
// for future use STL vector as container 
// typedef Vector<G4DynamicParticle> G4ParticleVector;
// #elseifdef RWT

typedef RWTPtrOrderedVector<G4DynamicParticle> G4ParticleVector;

#endif
