// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicsTable.hh,v 1.4 2001-01-09 01:18:50 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
//
//      Modified 01 March 1996, K. Amako
// ------------------------------------------------------------

#ifndef G4PhysicsTable_h
#define G4PhysicsTable_h 1

#include "g4rw/tpordvec.h"
#include "globals.hh"
class     G4PhysicsVector;

typedef G4RWTPtrOrderedVector<G4PhysicsVector> G4PhysicsTable;

#include "G4PhysicsVector.hh"
#endif


