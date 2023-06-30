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
// Declare custom type for fortran <-> C/C++ conversion.
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA
#ifndef TYPES_H
#define TYPES_H


#ifdef __cplusplus
extern "C" {
#else
// At least C99
#include <stdint.h>
#endif

// ***************************************************************************
// BOOLEAN
// ***************************************************************************
// Fortran LOGICAL typically takes 4 bytes by default,
// while C99+ / C++ bolean typically takes 1 byte.
//
// Hence, one could directly convert any Fortran LOGICAL to an int32_t,
// but for obvious readability and maintanability reasons
// (keep track of booleans),
// we declare here our custom 'logical' type.
//
// Important NB: Careful! The fact that a FORTRAN LOGICAL takes 4 bytes by default,
// is obviously NOT guaranteed by any FORTRAN standard though!
// A more portable approach would be to use ISO_C_BINDING 
// (but requires changes in FLUKA Fortran sources).

typedef int32_t logical;
// typedef instead of using, to be compatible with C.

#ifdef __cplusplus
}
#endif


#endif
#endif // G4_USE_FLUKA
