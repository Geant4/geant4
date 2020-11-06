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
// G4Types
//
// Definition of global GEANT4 native types

// Author: G.Cosmo (CERN), 1995
// --------------------------------------------------------------------
#ifndef G4TYPES_HH
#define G4TYPES_HH 1

// Fundamental definitions
#ifndef G4GMAKE
#  include "G4GlobalConfig.hh"
#endif

#ifdef WIN32
// Disable warning C4786 on WIN32 architectures:
// identifier was truncated to '255' characters
// in the debug information
//
#  pragma warning(disable : 4786)
//
// Define DLL export macro for WIN32 systems for
// importing/exporting external symbols to DLLs
//
#  if defined G4LIB_BUILD_DLL && !defined G4MULTITHREADED
#    define G4DLLEXPORT __declspec(dllexport)
#    define G4DLLIMPORT __declspec(dllimport)
#  else
#    define G4DLLEXPORT
#    define G4DLLIMPORT
#  endif
//
// Unique identifier for global module
//
#  if defined G4GLOB_ALLOC_EXPORT
#    define G4GLOB_DLL G4DLLEXPORT
#    define G4MTGLOB_DLL __declspec(dllexport)
#  else
#    define G4GLOB_DLL G4DLLIMPORT
#    define G4MTGLOB_DLL __declspec(dllimport)
#  endif
#else
#  define G4DLLEXPORT
#  define G4DLLIMPORT
#  define G4GLOB_DLL
#  define G4MTGLOB_DLL
#endif

#include <complex>

// Definitions for Thread Local Storage
//
#include "tls.hh"

// Typedefs to decouple from library classes
// Typedefs for numeric types
//
using G4double  = double;
using G4float   = float;
using G4int     = int;
using G4bool    = bool;
using G4long    = long;
using G4complex = std::complex<G4double>;

// Forward declation of void type argument for usage in direct object
// persistency to define fake default constructors
//
class __void__;

#endif /* G4TYPES_HH */
