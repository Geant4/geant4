// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VoxelLimits.hh,v 1.1 1999-01-07 16:07:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VoxelLimits
//
// Represents limitation/restrictions of space , where restrictions
// are only made perpendicular to the cartesian axes.
//
//
// Member functions:
//
// G4VoxelLimits()
//   Construct, with volume unrestricted
// ~G4VoxelLimits()
//   No actions.
// AddLimit(const EAxis pAxis, const G4double pMin,const G4double pMax)
//   Restict the volume to between specified min and max along the given axis.
//   Cartesian axes only, pMin<=pMax.
// G4double GetMaxXExtent() const
//   Return maximum x extent
// G4double GetMaxYExtent() const
//   Return maximum y extent
// G4double GetMaxZExtent() const
//   Return maximum z extent
// G4double GetMinXExtent() const
//   Return minimum x extent
// G4double GetMinYExtent() const
//   Return minimum y extent
// G4double GetMinZExtent() const
//   Return minimum z extent
// G4double GetMaxExtent(const EAxis pAxis) const
//   Return maximum extent of volume along specified axis.
// G4double GetMinExtent(const EAxis pAxis) const
//   Return maximum extent of volume along specified axis.
// G4bool IsLimited() const
//   Return true if limited along any axis
// G4bool IsLimited(const EAxis pAxis) const
//   Return true if the specified axis is resticted/limited.
// G4bool IsXLimited() const
//   Return true if the x axis is limited
// G4bool IsYLimited() const
//   Return true if the y axis is limited
// G4bool IsZLimited() const
//   Return true if the z axis is limited
//
// G4bool ClipToLimits(G4ThreeVector& pStart,G4ThreeVector& pEnd)
//   Clip the line segment pStart->pEnd to the volume described by the
//   current limits. Return true if the line remains after clipping,
//   else false, and leave the vectors in an undefined state.
//
// G4bool Inside(const G4ThreeVector& pVec) const
//   Return true if the specified vector is inside/on boundaries
//   of limits
//
// G4int OutCode(const G4ThreeVector& pVec) const
//   Calculate the `outcode' for the specified vector.
//   Intended for use during clipping against the limits
//   The bits are set given following conditions:
//     0      pVec.x()<fxAxisMin && IsXLimited()
//     1      pVec.x()>fxAxisMax && IsXLimited()
//     2      pVec.y()<fyAxisMin && IsYLimited()
//     3      pVec.y()>fyAxisMax && IsYLimited()
//     4      pVec.z()<fzAxisMin && IsZLimited()
//     5      pVec.z()>fzAxisMax && IsZLimited()
//
// Member data:
//
// G4double fxAxisMin,fxAxisMax
// G4double fyAxisMin,fyAxisMax
// G4double fzAxisMin,fzAxisMax
//   The min and max values along each axis. +-kInfinity if not restricted
//
//
// operators:
//
// ostream& operator << (ostream& os, const G4VoxelLimits& pLim);
//
// Print the limits to the stream in the form:
//  "{(xmin,xmax) (ymin,ymax) (zmin,zmax)}" Replace (xmin,xmax) by (-,-)
//  when not limited.
//
// Notes:
//
// Beware no break statements after returns in switch(pAxis)s
//
// History:
// 13.07.95 P.Kent Initial version.

#ifndef G4VOXELLIMITS_HH
#define G4VOXELLIMITS_HH

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"

#include <assert.h>

class ostream;

class G4VoxelLimits
{
public:
// Constructor - initialise to be unlimited
    G4VoxelLimits() : fxAxisMin(-kInfinity),fxAxisMax(kInfinity),
                      fyAxisMin(-kInfinity),fyAxisMax(kInfinity),
                      fzAxisMin(-kInfinity),fzAxisMax(kInfinity)
    {;}

//    G4VoxelLimits(const G4VoxelLimits& v);

// Destructor
    ~G4VoxelLimits() {;}

// Further restict limits
    void AddLimit(const EAxis pAxis, const G4double pMin,const G4double pMax);

// Return appropriate max limit
    G4double GetMaxXExtent() const
    {
	return fxAxisMax;
    }
    G4double GetMaxYExtent() const
    {
	return fyAxisMax;
    }
    G4double GetMaxZExtent() const
    {
	return fzAxisMax;
    }

// Return appropriate min limit
    G4double GetMinXExtent() const
    {
	return fxAxisMin;
    }
    G4double GetMinYExtent() const
    {
	return fyAxisMin;
    }
    G4double GetMinZExtent() const
    {
	return fzAxisMin;
    }

// Return specified max limit
    G4double GetMaxExtent(const EAxis pAxis) const
    {
	if (pAxis==kXAxis)
	    {
		return GetMaxXExtent();
	    }
	else if (pAxis==kYAxis)
	    {
		return GetMaxYExtent();
	    }
	else 
	    {
		assert(pAxis==kZAxis);
		return GetMaxZExtent();
	    }
    }

//Return min limit
    G4double GetMinExtent(const EAxis pAxis) const
    {
	if (pAxis==kXAxis)
	    {
		return GetMinXExtent();
	    }
	else if (pAxis==kYAxis)
	    {
		return GetMinYExtent();
	    }
	else 
	    {
		assert(pAxis==kZAxis);
		return GetMinZExtent();
	    }
    }

// Return true if x axis is limited
    G4bool IsXLimited() const
    {
	return (fxAxisMin==-kInfinity&&fxAxisMax==kInfinity) ? false : true;
    }
// Return true if y axis is limited
    G4bool IsYLimited() const
    {
	return (fyAxisMin==-kInfinity&&fyAxisMax==kInfinity) ? false : true;
    }
// Return true if z axis is limited
    G4bool IsZLimited() const
    {
	return (fzAxisMin==-kInfinity&&fzAxisMax==kInfinity) ? false : true;
    }

// Return true if limited along any axis
    G4bool IsLimited() const
    {
	return (IsXLimited()||IsYLimited()||IsZLimited());
    }

// Return true if specified axis is limited
    G4bool IsLimited(const EAxis pAxis) const
    {
	if (pAxis==kXAxis)
	    {
		return IsXLimited();
	    }
	else if (pAxis==kYAxis)
	    {
		return IsYLimited();
	    }
	else 
	    {
		assert(pAxis==kZAxis);
		return IsZLimited();
	    }
    }

    G4bool ClipToLimits(G4ThreeVector& pStart,G4ThreeVector& pEnd) const;

// Return true if specified vector is inside/on boundaries of limits
    G4bool Inside(const G4ThreeVector& pVec) const
    {
	return ((GetMinXExtent()<=pVec.x()) &&
		(GetMaxXExtent()>=pVec.x()) &&
		(GetMinYExtent()<=pVec.y()) &&
		(GetMaxYExtent()>=pVec.y()) &&
		(GetMinZExtent()<=pVec.z()) &&
		(GetMaxZExtent()>=pVec.z()) ) ? true : false;
    }

    G4int OutCode(const G4ThreeVector& pVec) const;

private:
    G4double fxAxisMin,fxAxisMax;
    G4double fyAxisMin,fyAxisMax;
    G4double fzAxisMin,fzAxisMax;
};

ostream& operator << (ostream& os, const G4VoxelLimits& pLim);

#endif
