// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ios.hh,v 1.3 1999-11-16 17:40:51 gcosmo Exp $
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

#if defined(OO_DDL_TRANSLATION)
/*
 * stdlib needs to be included before iostream.h
 * during oddlx runs to work around a parser
 * problem of AIX ooddlx v4.0.2 with some versions of
 * AIX system header files.
 */
#include <stdlib.h>
#endif

#include <iostream.h>

  extern ostream G4cout;
  extern ostream G4cerr;

#endif
