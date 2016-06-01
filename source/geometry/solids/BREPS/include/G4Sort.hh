// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sort.hh,v 2.1 1998/10/20 16:31:33 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
//
//  File: G4Sort.h
//  Author:  Alan Breakstone
//
//  Description
//
//	Routines to G4Sort arG4Rays of various kinds of numbers
//

#ifndef __SORT_H
#define __SORT_H
#include "globals.hh"

void G4Sort_double( G4double [], int, int );

void swap_double(   G4double [], int, int );

void G4Sort_float( float [], int, int );

void swap_float(   float [], int, int );

#endif
