// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ios.hh,v 1.5 1999-12-15 18:05:18 gracia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// G4ios.hh
//
// ---------------------------------------------------------------
#ifndef included_G4ios
#define included_G4ios

#include "G4Types.hh"

#if defined(OO_DDL_TRANSLATION)
/*
 * stdlib needs to be included before iostream.h
 * during oddlx runs to work around a parser
 * problem of AIX ooddlx v4.0.2 with some versions of
 * AIX system header files.
 */
#include <stdlib.h>
#endif

#include "g4std/iostream"

extern G4std::ostream G4cout;
extern G4std::ostream G4cerr;

#define G4cin G4std::cin
#define G4endl G4std::endl

#endif
