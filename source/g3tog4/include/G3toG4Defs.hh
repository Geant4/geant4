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
#ifndef G3TOG4_DEFS_H
#define G3TOG4_DEFS_H

#ifdef WIN32
  //
  // Define DLL export macro for WIN32 systems for
  // importing/exporting external symbols to DLLs
  //
  #if defined G4LIB_BUILD_DLL
    #if defined G3TOG4_EXPORT
      #define G3G4DLL_API __declspec( dllexport )
    #else
      #define G3G4DLL_API __declspec( dllimport )
    #endif
  #else
    #define G3G4DLL_API
  #endif
#else
  #define G3G4DLL_API
#endif

#endif /* G3TOG4_DEFS_H */
