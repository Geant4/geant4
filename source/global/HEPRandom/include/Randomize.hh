// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Randomize.hh,v 2.1 1998/09/23 23:11:15 gcosmo Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef randomize_h
#define randomize_h 1

#include <CLHEP/Random/Randomize.h>

#define G4UniformRand() HepRandom::getTheEngine()->flat()

#endif // randomize_h 
