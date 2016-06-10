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
//
// $Id: $
//
//
// Defines for Windows DLLs import/export
//

#ifndef TRKGDEFS_HH
#define TRKGDEFS_HH

#include "G4Types.hh"

#ifdef WIN32
  //
  // Unique identifier for global module
  //
  #if defined G4TRACKING_ALLOC_EXPORT
    #define G4TRACKING_DLL G4DLLEXPORT
  #else
    #define G4TRACKING_DLL G4DLLIMPORT
  #endif
#else
  #define G4TRACKING_DLL
#endif

#endif /* G4TRKGDEFS_HH */
