// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sort.hh,v 1.3 2000-08-28 08:57:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Global functions for sorting and swaping arrays of various
// kinds of numbers.

// Author:  Alan Breakstone
// ----------------------------------------------------------------------
#ifndef __SORT_H
#define __SORT_H

#include "globals.hh"

void sort_double( G4double [], G4int, G4int );

void swap_double( G4double [], G4int, G4int );

#endif
