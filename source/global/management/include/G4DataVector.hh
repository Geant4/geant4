// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DataVector.hh,v 1.2 1999-11-11 10:47:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
// ------------------------------------------------------------

#ifndef G4DataVector_h
#define G4DataVector_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "g4rw/tvordvec.h"

typedef G4RWTValOrderedVector<G4double> G4DataVector;

#endif
