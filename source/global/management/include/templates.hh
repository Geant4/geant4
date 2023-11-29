//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Module defining platform dependent features and some useful utilities.

// Author: Gabriele Cosmo, 5 September 1995 - Created
// --------------------------------------------------------------------
#ifndef templates_hh
#define templates_hh 1

#include <climits>
#include <limits>

// If HIGH_PRECISION is defined to TRUE (ie. != 0) then the type "Float"
// is typedefed to "double". If it is FALSE (ie. 0) it is typedefed
// to "float".
//
#ifndef HIGH_PRECISION
#  define HIGH_PRECISION 1
#endif

#if HIGH_PRECISION
using Float = double;
#else
using Float = float;
#endif

// Following values have been taken from limits.h
// and temporarly defined for portability on HP-UX.

#ifndef DBL_MIN /* Min decimal value of a double */
#  define DBL_MIN std::numeric_limits<double>::min()  // 2.2250738585072014e-308
#endif

#ifndef DBL_DIG /* Digits of precision of a double */
#  define DBL_DIG std::numeric_limits<double>::digits10  // 15
#endif

#ifndef DBL_MAX /* Max decimal value of a double */
#  define DBL_MAX std::numeric_limits<double>::max()  // 1.7976931348623157e+308
#endif

#ifndef DBL_EPSILON
#  define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif  // 2.2204460492503131e-16

#ifndef FLT_MIN /* Min decimal value of a float */
#  define FLT_MIN std::numeric_limits<float>::min()  // 1.17549435e-38F
#endif

#ifndef FLT_DIG /* Digits of precision of a float */
#  define FLT_DIG std::numeric_limits<float>::digits10  // 6
#endif

#ifndef FLT_MAX /* Max decimal value of a float */
#  define FLT_MAX std::numeric_limits<float>::max()  // 3.40282347e+38F
#endif

#ifndef FLT_EPSILON
#  define FLT_EPSILON std::numeric_limits<float>::epsilon()
#endif  // 1.192092896e-07F

#ifndef MAXFLOAT /* Max decimal value of a float */
#  define MAXFLOAT std::numeric_limits<float>::max()  // 3.40282347e+38F
#endif

#ifndef INT_MAX /* Max decimal value of a int */
#  define INT_MAX std::numeric_limits<int>::max()  // 2147483647
#endif

#ifndef INT_MIN /* Min decimal value of a int */
#  define INT_MIN std::numeric_limits<int>::min()  // -2147483648
#endif

#ifndef LOG_EKIN_MIN /* Min value of the natural logarithm of kin. energy. */
#  define LOG_EKIN_MIN -30
#endif

//---------------------------------

template <class T>
inline void G4SwapPtr(T*& a, T*& b)
{
  T* tmp = a;
  a      = b;
  b      = tmp;
}

template <class T>
inline void G4SwapObj(T* a, T* b)
{
  T tmp = *a;
  *a    = *b;
  *b    = tmp;
}

//-----------------------------

#ifndef G4_SQR_DEFINED
#  define G4_SQR_DEFINED
#  ifdef sqr
#    undef sqr
#  endif

template <class T>
inline T sqr(const T& x)
{
  return x * x;
}
#endif

inline int G4lrint(double ad)
{
  return (int)std::lrint(ad);
}

//-----------------------------

//  Use the following function to get rid of "unused parameter" warnings
//  Example:
//
//      #ifdef SOME_CONDITION
//          void doSomething(int val)
//          {
//              something = val;
//          }
//      #else
//          void doSomething(int)
//          { }
//      #endif
//
//  can be simplified to:
//
//          void doSomething(int val)
//          {
//      #ifdef SOME_CONDITION
//              something = val;
//      #else
//              G4ConsumeParameters(val);
//      #endif
//          }
//
//  or:
//
//          void doSomething(int val)
//          {
//      #ifdef SOME_CONDITION
//              something = val;
//      #endif
//              // function call does nothing -- will be "optimized" out
//              G4ConsumeParameters(val);
//          }
//
template <typename... _Args>
inline void G4ConsumeParameters(_Args&&...)
{}

#endif  // templates_hh
