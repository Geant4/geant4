//
// G4SolidExtentList.hh
//
// Declaration of a list of (voxel) extents along one axis
//
// This utility class is designed for one specific purpose: to
// calculate the extent of a CSG solid for a voxel
// (G4VSolid::CalculateExtent). 
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef G4SolidExtentList_hh
#define G4SolidExtentList_hh

#include "globals.hh"

#include "G4ClippablePolygon.hh"

class G4SolidExtentList {
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

	G4ClippablePolygon minSurface,		// Minimum surface within limits
			   maxSurface,		// Maximum
			   minAbove,		// Minimum surface totally above max limit
			   maxBelow;		// Maximum surface totally below min limit
};


#endif
