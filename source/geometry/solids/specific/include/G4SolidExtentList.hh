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
// $Id: G4SolidExtentList.hh,v 1.4 2001-07-11 10:00:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//   to calculate the extent of a CSG solid for a voxel
//   (G4VSolid::CalculateExtent). 

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4SolidExtentList_hh
#define G4SolidExtentList_hh

#include "globals.hh"

#include "G4ClippablePolygon.hh"

class G4SolidExtentList
{
  public:
	
	G4SolidExtentList();
	G4SolidExtentList( const EAxis targetAxis, const G4VoxelLimits &voxelLimits );
	~G4SolidExtentList();


	void AddSurface( const G4ClippablePolygon &surface );

	G4bool GetExtent( G4double &min, G4double &max ) const;

	protected:
	
	EAxis	 axis;		// Target axis
	G4bool   limited;	// True if limited
	G4double minLimit;	// ... min limit
	G4double maxLimit;	// ... max limit

	G4ClippablePolygon minSurface,  // Minimum surface within limits
			   maxSurface,  // Maximum
			   minAbove,    // Minimum surface totally above max limit
			   maxBelow;    // Maximum surface totally below min limit
};


#endif
