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
// $Id: G4Types.hh,v 1.13 2005/11/09 10:26:23 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// GEANT4 native types
//
// History:

#ifndef G4TYPES_HH
#define G4TYPES_HH

#ifdef WIN32
  // Disable warning C4786 on WIN32 architectures:
  // identifier was truncated to '255' characters
  // in the debug information
  //
  #pragma warning ( disable : 4786 )
  //
  // Define DLL export macro for WIN32 systems for
  // importing/exporting external symbols to DLLs
  //
  #if defined G4LIB_BUILD_DLL
    #define G4DLLEXPORT __declspec( dllexport )
    #define G4DLLIMPORT __declspec( dllimport )
  #else
    #define G4DLLEXPORT
    #define G4DLLIMPORT
  #endif
#else
  #define G4DLLEXPORT
  #define G4DLLIMPORT
#endif

#include <complex>

// Typedefs to decouple from library classes
// Typedefs for numeric types
//
typedef double G4double;
typedef float G4float;
typedef int G4int;
typedef bool G4bool;
typedef long G4long;
typedef std::complex<G4double> G4complex;

// Forward declation of void type argument for usage in direct object
// persistency to define fake default constructors
//
class __void__;

#endif /* G4TYPES_HH */
