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
// $Id: voxeldefs.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
//
// Voxel Optimisation Constants

// History:
// 12.02.99 S.Giani made numerical values consistent with Geant3.21
// 13.08.95 P.Kent Created separate file
// --------------------------------------------------------------------
#ifndef VOXELDEFS_HH
#define VOXELDEFS_HH

#include "G4Types.hh"

// Hard limit on no. voxel nodes per given header
const G4int kMaxVoxelNodes=1000;  // PK chose 2000, Geant 3.21 used 1000

const G4int kMinVoxelVolumesLevel1=2; // Only begin to make voxels if >=
				      // this no of daughters
const G4int kMinVoxelVolumesLevel2=3; // Only make second level of refinement
				      // if >= this no of volumes in
                                      // 1st level node
const G4int kMinVoxelVolumesLevel3=4; // Only make third level of refinement
				      // if >= this no of volumes in
                                      // 2nd level node
#endif
