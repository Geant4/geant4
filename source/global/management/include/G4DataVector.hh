// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DataVector.hh,v 1.1 1999-01-07 16:09:00 gunter Exp $
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
#include <rw/tvordvec.h>

typedef RWTValOrderedVector<G4double> G4DataVector;

#endif
