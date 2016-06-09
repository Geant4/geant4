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
// $Id: G4Types.hh,v 1.6 2003/06/06 16:17:14 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
//
// GEANT4 native types
//
// History:

#ifndef G4TYPES_HH
#define G4TYPES_HH

// Disable warning C4786 on WIN32 architectures: identifier was truncated
// to '255' characters in the debug information
//
#ifdef WIN32
  #pragma warning ( disable : 4786 )
#endif

// Disable deprecated warnings for usage of strstream on Linux
// architectures with gcc >= 3.0 release
//
#if (__GNUC__==3) && (__GNUC_MINOR__>0)
  #undef __DEPRECATED
#endif

#include <CLHEP/config/CLHEP.h>
#include <complex>

// Define G4std namespace for standard std.
// Now obsolete. Will be removed from release 6.0.
//
#define G4std std

// Typedefs to decouple from library classes
// Typedefs for numeric types
//
typedef double G4double;
typedef float G4float;
typedef int G4int;
typedef bool G4bool;
typedef long G4long;
typedef std::complex<G4double> G4complex;

#endif /* G4TYPES_HH */
