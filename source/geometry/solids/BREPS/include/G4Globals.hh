// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Globals.hh,v 1.1 1999-01-07 16:07:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*  /usr/local/gismo/repo/support/globals.h,v 1.8 1994/04/18 18:29:03 atwood Exp  */
//  File: globals.h
//  Author: Alan Breakstone
//
//  Description
//
//      Defines global variables and declarations for Gismo
//
//

#ifndef __GLOBALS_H
#define __GLOBALS_H

//
//  Define a C preprocessor constant for a Scale factor to apply to 
//  various dimensionless tests in the geometry routines which test
//  if a point is on a surface.  This number is an effective thickness
//  of a surface divided by a relevant dimension, such as the radius of
//  a cylinder.  The default value is 0.0001.
#define SURFACE_PRECISION 0.0001 
//
//  Define a C preprocessor constant for the maximum number of turns
//  allowed for a Helix, which is used in some of the geometry routines
//  to limit the size of for or while loops.  The default value is 50.
#define HELIX_MAX_TURNS 50

//  Define some geometric constants of use
//  These should be gotten via math.h ... see M_PI etc....  
//#define PI 3.14159265358979323846
//#define TWO_PI 6.2831853071795862

#endif
