//
// G4SolidExtentList.hh
//
// Declaration of a list of (voxel) extents along one axis
//
// This utility class is designed for one specific purpose: to
// calculate the extent of a CSG solid for a voxel
// (G4VSolid::CalculateExtent). It consists of a sorted, single-linked
// list of objects which are characterized by two values: 
//         1. A double value (indicating the distance along
//            an axis associated with a solid surface
//         2. Whether this value is asociated with an "outgoing"
//            surface.
// Outgoing means that the normal of the surface (pointing outwards
// from the solid) is aligned with the axis (dot product with axis
// is positive).
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

class G4SolidExtentList {
	public:
	
	G4SolidExtentList();
	~G4SolidExtentList();


	void AddSurface( const G4double aMin, const G4double aMax,
			 const G4double aMinNormal, const G4double aMaxNormal,
			 const G4double tolerance );

	G4bool GetExtent( G4double &min, G4double &max ) const;
	G4bool UpdateLimitedExtent( G4double &min, G4double &max ) const;

	protected:

	G4double min, max;		// The current min and max
	G4double minNormal,
		 maxNormal;		// ...and the corresponding normals
};


#endif
