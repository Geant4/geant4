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
// $Id: G4SolidExtentList.cc,v 1.3 2002-10-28 11:47:53 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


//
// Constructor (default)
//
G4SolidExtentList::G4SolidExtentList() 
{
  axis = kZAxis;
  limited = false;
  minLimit = -DBL_MAX;
  maxLimit = +DBL_MAX;
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
    minLimit = -DBL_MAX;
    maxLimit = +DBL_MAX;
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

