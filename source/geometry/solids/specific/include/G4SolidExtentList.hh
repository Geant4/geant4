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
// $Id: G4SolidExtentList.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4SolidExtentList
//
// Class description:
//
//   Defines a list of (voxel) extents along one axis.
//
//   This utility class is designed for one specific purpose:
//   to calculate the extent of a CSG-like solid for a voxel
//   (G4VSolid::CalculateExtent). 

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4SolidExtentList_hh
#define G4SolidExtentList_hh

#include "G4Types.hh"

#include "G4ClippablePolygon.hh"

class G4SolidExtentList
{
  public:
  
  G4SolidExtentList();
  G4SolidExtentList( const EAxis targetAxis,
                     const G4VoxelLimits &voxelLimits );
  ~G4SolidExtentList();


  void AddSurface( const G4ClippablePolygon &surface );

  G4bool GetExtent( G4double &min, G4double &max ) const;

  protected:
  
  EAxis    axis;     // Target axis
  G4bool   limited;  // True if limited
  G4double minLimit; // ... min limit
  G4double maxLimit; // ... max limit

  G4ClippablePolygon minSurface,  // Minimum surface within limits
                     maxSurface,  // Maximum
                     minAbove,    // Minimum surface totally above max limit
                     maxBelow;    // Maximum surface totally below min limit
};

#endif
