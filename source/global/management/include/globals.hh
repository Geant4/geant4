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
// $Id: globals.hh,v 1.26 2005/11/04 08:18:51 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
// 15.12.99 G.Garcia - Included min, max definitions for NT with ISO standard
// 15.06.01 G.Cosmo - Removed cbrt() definition

#ifndef GLOBALS_HH
#define GLOBALS_HH

#include "G4ios.hh"

#ifndef FALSE
  #define FALSE 0
#endif
#ifndef TRUE
  #define TRUE 1
#endif

#include <algorithm>  // Retrieve definitions of min/max

// Include base types
#include "G4Types.hh"

// Get definition of G4String
#include "G4String.hh"

// Includes some additional definitions: sqr, G4SwapPtr, G4SwapObj.
#include "templates.hh"

// Includes Physical Constants and System of Units
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Global error function
#include "G4ExceptionSeverity.hh"
void G4Exception(const char* issure,
                 const char* errorCode,
                             G4ExceptionSeverity severity,
                 const char* comments);
void G4Exception(const char* s=0);
void G4Exception(std::string s);
void G4Exception(G4String s);

#endif /* GLOBALS_HH */

