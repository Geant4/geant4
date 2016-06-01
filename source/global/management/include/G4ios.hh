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

#ifdef G4STREAM
  extern ostream G4cout;
  extern ostream G4cerr;
#else
  #define G4cout cout
  #define G4cerr cerr
#endif

#endif
