// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: globals.hh,v 1.15 2000-01-19 11:12:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Global Constants and typedefs
//
// History:
// 30.06.95 P.Kent - Created
// 16.02.96 G.Cosmo - Added inclusion of "templates.hh"
// 03.03.96 M.Maire - Added inclusion of "G4PhysicalConstants.hh"
// 08.11.96 G.Cosmo - Added cbrt() definition and G4ApplicationState enum type
// 29.11.96 G.Cosmo - Added typedef of HepBoolean to G4bool
// 22.10.97 M.Maire - Moved PhysicalConstants at the end of the file
// 04.12.97 G.Cosmo,E.Tcherniaev - Migrated to CLHEP
// 26.08.98 J.Allison,E.Tcherniaev - Introduced min/max/sqr/abs functions
// 22.09.98 G.Cosmo - Removed min/max/sqr/abs functions and replaced with
//                    inclusion of CLHEP/config/TemplateFunctions.h for CLHEP-1.3
// 15.12.99 G.Gracia - Included min, max definitions for NT with ISO standard

#ifndef GLOBALS_HH
#define GLOBALS_HH

#include "G4ios.hh"

// Undefine possible existing min/max/sqr/abs macros first
// (temporary solution)
#ifdef min
  #undef min
#endif
#ifdef max
  #undef max
#endif
#ifdef sqr
  #undef sqr
#endif
#ifdef abs
  #undef abs
#endif

#include "g4std/algorithm"
#ifndef CLHEP_MAX_MIN_DEFINED
  #define CLHEP_MAX_MIN_DEFINED
#endif

#if defined(WIN32) && defined(G4USE_STD_NAMESPACE)
// For NT with Native STL (used in ISO standard mode)
// templated functions min and max should be _MIN _MAX
  #define min _MIN
  #define max _MAX
#endif

// min, max, abs and sqr are in TemplateFunctions.h.
// Includes also CLHEP.h with typedef for numeric types and
// implicit inclusions of <stdlib.h>, <limits.h>, <math.h>.
#include <CLHEP/config/TemplateFunctions.h>

// Include base types
#include "G4Types.hh"

// Get definition of G4String
#include "g4rw/cstring.h"

// Includes some additional definitions
#include "templates.hh"

// cbrt() function - define G4_NO_CBRT if the function is not available
#ifdef G4_NO_CBRT
  static double cbrt(double x) { return pow(x,1./3.); }
#endif

// System of Units and Physical Constants
////#include <CLHEP/Units/PhysicalConstants.h>
#include "PhysicalConstants.h"

// Global error function
void G4Exception(const char* s=0);
void G4Exception(G4std::string s);
void G4Exception(G4String s);

#endif /* GLOBALS_HH */

