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
// $Id: G4SolidExtentList.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4SolidExtentList.cc
//
// Implementation of a list of (voxel) extents along one axis
//
// --------------------------------------------------------------------

#include "G4SolidExtentList.hh"
#include "G4VoxelLimits.hh"
#include "G4GeometryTolerance.hh"


//
// Constructor (default)
//
G4SolidExtentList::G4SolidExtentList() 
{
  axis = kZAxis;
  limited = false;
  minLimit = -INT_MAX/2;
  maxLimit =  INT_MAX/2;
}


//
// Constructor (limited case)
//
G4SolidExtentList::G4SolidExtentList( const EAxis targetAxis,
                                      const G4VoxelLimits &voxelLimits )
{
  axis = targetAxis;
  
  limited = voxelLimits.IsLimited( axis );
  if (limited)
  {
    minLimit = voxelLimits.GetMinExtent( axis );
    maxLimit = voxelLimits.GetMaxExtent( axis );
  }
  else
  {
    minLimit = -INT_MAX/2;
    maxLimit =  INT_MAX/2;
  }
}


//
// Destructor
//
G4SolidExtentList::~G4SolidExtentList()
{
}


//
// AddSurface
//
//
void G4SolidExtentList::AddSurface( const G4ClippablePolygon &surface )
{
  //
  // Keep track of four surfaces
  //
  G4double min, max;
  
  surface.GetExtent( axis, min, max );
  
  if (min > maxLimit)
  {
    //
    // Nearest surface beyond maximum limit
    //
    if (surface.InFrontOf(minAbove,axis)) minAbove = surface;
  }
  else if (max < minLimit)
  {
    //
    // Nearest surface below minimum limit
    //
    if (surface.BehindOf(maxBelow,axis)) maxBelow = surface;
  }
  else
  {
    //
    // Max and min surfaces inside
    //
    if (surface.BehindOf(maxSurface,axis)) maxSurface = surface;
    if (surface.InFrontOf(minSurface,axis)) minSurface = surface;
  }
}



//
// GetExtent
//
// Return extent after processing all surfaces
//
G4bool G4SolidExtentList::GetExtent( G4double &min, G4double &max ) const
{
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()
                           ->GetSurfaceTolerance();
  //
  // Did we have any surfaces within the limits?
  //
  if (minSurface.Empty())
  {
    //
    // Nothing! Do we have anything above?
    //
    if (minAbove.Empty()) return false;
    
    //
    // Yup. Is it facing inwards?
    //
    if (minAbove.GetNormal().operator()(axis) < 0) return false;
    
    //
    // No. We must be entirely within the solid
    //
    max = maxLimit + kCarTolerance;
    min = minLimit - kCarTolerance;
    return true;
  }
  
  //
  // Check max surface
  //
  if (maxSurface.GetNormal().operator()(axis) < 0)
  {
    //
    // Inward facing: max limit must be embedded within solid
    //
    max = maxLimit + kCarTolerance;
  }
  else
  {
    G4double sMin, sMax;
    maxSurface.GetExtent( axis, sMin, sMax );
    max = ( (sMax > maxLimit) ? maxLimit : sMax ) + kCarTolerance;
  }
  
  //
  // Check min surface
  //
  if (minSurface.GetNormal().operator()(axis) > 0)
  {
    //
    // Inward facing: max limit must be embedded within solid
    //
    min = minLimit - kCarTolerance;
  }
  else
  {
    G4double sMin, sMax;
    minSurface.GetExtent( axis, sMin, sMax );
    min = ( (sMin < minLimit) ? minLimit : sMin ) - kCarTolerance;
  }
  
  return true;
}

