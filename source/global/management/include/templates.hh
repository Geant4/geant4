// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: templates.hh,v 1.4 1999-12-15 14:50:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// -*- C++ -*-
//
// -----------------------------------------------------------------------
// This file should define some platform dependent features and some
// useful utilities.
// -----------------------------------------------------------------------

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
// Gabriele Cosmo - Minor change: 08/02/1996
// Gabriele Cosmo - Added DBL_MIN, FLT_MIN, DBL_DIG,
//                        DBL_MAX, FLT_DIG, FLT_MAX  : 12/04/1996
// Gabriele Cosmo - Removed boolean enum definition : 29/11/1996
// Gunter Folger  - Added G4SwapPtr() and G4SwapObj() : 31/07/1997
// Gabriele Cosmo - Adapted signatures of min(), G4std::max() to
//                  STL's ones, thanks to E.Tcherniaev : 31/07/1997
// Gabriele Cosmo,
// Evgueni Tcherniaev - Migrated to CLHEP: 04/12/1997 
// =======================================================================

#ifndef templates_h
#define templates_h 1

//
// If HIGH_PRECISION is defined to TRUE (ie. != 0) then the type "Float"
// is typedefed to "double". If it is FALSE (ie. 0) it is typedefed
// to "float".
//
#ifndef HIGH_PRECISION
#define HIGH_PRECISION 1
#endif

#if HIGH_PRECISION
typedef double Float;
#else
typedef float Float;
#endif

// Following values have been taken from limits.h
// and temporarly defined for portability on HP-UX.

#ifndef DBL_MIN   /* Min decimal value of a double */
#define DBL_MIN   2.2250738585072014e-308
#endif

#ifndef DBL_DIG   /* Digits of precision of a double */
#define DBL_DIG   15
#endif

#ifndef DBL_MAX   /* Max decimal value of a double */
#define DBL_MAX   1.7976931348623157e+308
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON   2.2204460492503131e-16
#endif

#ifndef FLT_MIN   /* Min decimal value of a float */
#define FLT_MIN   1.17549435e-38F
#endif

#ifndef FLT_DIG   /* Digits of precision of a float */
#define FLT_DIG   6
#endif

#ifndef FLT_MAX   /* Max decimal value of a float */
#define FLT_MAX   3.40282347e+38F
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON   1.192092896e-07F
#endif

#ifndef MAXFLOAT   /* Max decimal value of a float */
#define MAXFLOAT   3.40282347e+38F
#endif

//---------------------------------

template <class T>
inline void G4SwapPtr(T* a, T* b) {
  T* tmp=a;
  a = b;
  b = tmp;
}

template <class T>
inline void G4SwapObj(T* a, T* b) {
  T tmp= *a;
  *a = *b;
  *b = tmp;
}

//-----------------------------

#endif // templates_h
