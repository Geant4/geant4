// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Randomize.hh,v 1.1 1999-01-07 16:08:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef randomize_h
#define randomize_h 1

#include <CLHEP/Random/Randomize.h>

#define G4UniformRand() HepRandom::getTheEngine()->flat()

#endif // randomize_h 
