// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Types.hh,v 1.1 1999-05-24 18:23:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 native types
//
// History:

#ifndef G4TYPES_HH
#define G4TYPES_HH

#include <CLHEP/config/CLHEP.h>

// Typedefs to decouple from library classes
// Typedefs for numeric types
// [NOTE: Will in future need to be made more sophisticated]
typedef HepDouble G4double;
typedef HepFloat G4float;
typedef HepInt G4int;
#ifdef G4_HAVE_BOOL
  typedef bool G4bool;
#else
  typedef HepBoolean G4bool;
#endif
typedef long G4long;

#endif /* G4TYPES_HH */

