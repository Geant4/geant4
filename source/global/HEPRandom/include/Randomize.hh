// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Randomize.hh,v 1.2 1999-11-16 17:31:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef randomize_h
#define randomize_h 1

#include <CLHEP/Random/Randomize.h>

#define G4UniformRand() HepRandom::getTheEngine()->flat()

#endif // randomize_h 
