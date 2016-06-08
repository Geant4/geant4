// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4o.h,v 1.3 1999/12/15 14:48:41 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
/* +---------------------- Copyright notice -------------------------------+ */
/* | Copyright (C) 1995, Guy Barrand, LAL Orsay, (barrand@lal.in2p3.fr)    | */
/* |   Permission to use, copy, modify, and distribute this software       | */
/* |   and its documentation for any purpose and without fee is hereby     | */
/* |   granted, provided that the above copyright notice appear in all     | */
/* |   copies and that both that copyright notice and this permission      | */
/* |   notice appear in supporting documentation.  This software is        | */
/* |   provided "as is" without express or implied warranty.               | */
/* +---------------------- Copyright notice -------------------------------+ */
#ifndef G4o_h
#define G4o_h
 
#include <OShell.h>
#ifdef __cplusplus
extern "C"{
#endif
void   G4oAddCommands             (void*);
void   G4oExecuteScript           (char*);
OShell G4oGetShell                ();
#ifdef __cplusplus
}
#endif

#endif /*G4o_h*/
