// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sort.hh,v 1.1 1999-01-07 16:07:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
