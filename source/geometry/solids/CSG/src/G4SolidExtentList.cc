//
// G4SolidExtentList.cc
//
// Implementation of a list of (voxel) extents along one axis
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4SolidExtentList.hh"


//
// Constructor
//
G4SolidExtentList::G4SolidExtentList() 
{
	max = -DBL_MAX;
	min = +DBL_MAX;
}


//
// Destructor
//
G4SolidExtentList::~G4SolidExtentList() {;}



//
// AddExtent
//
// Add a new extent. Arguments:
//	aMin	- The minimum of the clipped surface along the target axis
//	aMax	- The maximum of the clipped surface along the target axis
//	aNinNormal	- The dot product of the target axis and the 
//			  normal of the surface at the minimum point
//	aNaxNormal	- The dot product of the target axis and the 
//			  normal of the surface at the maximum point
//	tolerance	- Surface tolerance
//
// Notes:
//	aMinNormal and aMaxNormal will be equal for a planer surface
//
void G4SolidExtentList::AddSurface( const G4double aMin, const G4double aMax, 
			            const G4double aMinNormal, const G4double aMaxNormal,
				    const G4double tolerance )
{
	//
	// Decide if we have a new maximum or minimum, and
	// update the 
	//
	// Be very careful how you decide: we want to avoid updating
	// maxNormal if (due to roundoff problems)
	// a maximum belonging to a surface with a 
	// positive normal accidently falls below a maximum belonging
	// to a surface with a negative normal.
	//
	// We do, however, want to keep track of the absolute maximum
	// and minimum, to be as conservative as possible.
	//
	if (aMinNormal >= 0) {
		if (aMax > max-tolerance) maxNormal = aMaxNormal;
	}
	else {
		if (aMax > max+tolerance) maxNormal = aMaxNormal;
	}
	if (aMax > max) max = aMax;


	if (aMaxNormal <= 0) {
		if (aMin < min+tolerance) minNormal = aMinNormal;
	}
	else {
		if (aMin > min-tolerance) minNormal = aMinNormal;
	}
	if (aMin < min) min = aMin;
}



//
// GetExtent
//
// Return extent for the unlimited case
//
G4bool G4SolidExtentList::GetExtent( G4double &theMin, G4double &theMax ) const
{
	if (min > max) return false;
	
	theMin = min;
	theMax = max;
	
	return true;
}


//
// UpdateLimitedExtent
//
// Return extent limited by input values of min and max. If nothing
// is left of the solid within the input values, return false.
//
G4bool G4SolidExtentList::UpdateLimitedExtent( G4double &updatedMin, G4double &updatedMax ) const
{
	if (min > max) return false;
	
	//
	// Update max, but only if the normal of the surface that produced
	// the maximum is positive
	//
	if (maxNormal >= 0) updatedMax = max;
	
	//
	// Similar for minimum
	//
	if (minNormal <= 0) updatedMin = min;
	
	return true;
}

	
