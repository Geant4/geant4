//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: globals.hh,v 1.4 2001-07-11 09:59:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Global Constants and typedefs
//
// History:
// 30.06.95 P.Kent

#ifndef GLOBALS_HH
#define GLOBALS_HH

// Typedefs for numeric types
// [NOTE: Will in future need to be made more sophisticated]
typedef double G4double;
typedef float G4float;
typedef long G4long;
typedef int G4int;

// Typedefs to decouple from library classes
//#include <rw/cstring.h>
//typedef G4String G4String;

// Boolean - define G4_HAVE_BOOL if bool type available
#ifdef G4_HAVE_BOOL
typedef bool G4bool;
#else
typedef int G4bool;
const int false = 0;
const int true = 1;
//enum G4bool {false = 0, true = 1};
#endif

// Global error function
//void G4Exception(const char* s=0);

#endif





