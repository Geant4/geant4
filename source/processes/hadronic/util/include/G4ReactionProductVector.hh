// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReactionProductVector.hh,v 1.2 1999-11-11 15:37:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
//	History: first implementation, alternative to G4FastVector
//               less fast, but it has variable array length and checks boundaries
//	26th September, Chr. Voelcker
// ------------------------------------------------------------

#ifndef G4ReactionProductVector_h
#define G4ReactionProductVector_h 1

#include "globals.hh"
#include "G4ios.hh"

class G4ReactionProduct;
#include "g4rw/tpordvec.h"

// #ifdef STL
// //in future use STL vector as container of reaction products ...
// typedef Vector<G4ReactionProduct> G4ReactionProductVector;
// #elseifdef RWT

typedef G4RWTPtrOrderedVector<G4ReactionProduct> G4ReactionProductVector;

// #endif

#endif
