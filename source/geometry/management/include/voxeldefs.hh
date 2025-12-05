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
// Voxel Optimisation Constants

// Author: Paul Kent (CERN), 13.08.1995
// --------------------------------------------------------------------
#ifndef VOXELDEFS_HH
#define VOXELDEFS_HH

#include "G4Types.hh"

/** Hard limit on the number of voxel nodes per given header. */
const G4int kMaxVoxelNodes = 1000;  // PK chose 2000, Geant 3.21 used 1000

/** Only begin to make voxels if greater or equal to this number of daughters. */
const G4int kMinVoxelVolumesLevel1 = 2;

/** Only make second level of refinement if greater or equal to this number of
    volumes in the first level node. */
const G4int kMinVoxelVolumesLevel2 = 3;

/** Only make third level of refinement if greater or equal to this number of
    volumes in the second level node. */
const G4int kMinVoxelVolumesLevel3 = 4;

#endif
