// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VoxelLimits.hh,v 1.4 2000-04-20 16:49:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VoxelLimits
//
// Class description:
//
// Represents limitation/restrictions of space, where restrictions
// are only made perpendicular to the cartesian axes.
//
//
// Member data:
//
// G4double fxAxisMin,fxAxisMax
// G4double fyAxisMin,fyAxisMax
// G4double fzAxisMin,fzAxisMax
//   - The min and max values along each axis. +-kInfinity if not restricted.
//
//
// Notes:
//
// Beware no break statements after returns in switch(pAxis)s.

// History:
// 13.07.95 P.Kent Initial version.

#ifndef G4VOXELLIMITS_HH
#define G4VOXELLIMITS_HH

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "g4std/iostream"

#include <assert.h>

class G4VoxelLimits
{
  public: // with description
  
    G4VoxelLimits() : fxAxisMin(-kInfinity),fxAxisMax(kInfinity),
                      fyAxisMin(-kInfinity),fyAxisMax(kInfinity),
                      fzAxisMin(-kInfinity),fzAxisMax(kInfinity)  {;}
      // Constructor - initialise to be unlimited. Volume unrestricted.

    ~G4VoxelLimits() {;}
      // Destructor. No actions.

    void AddLimit(const EAxis pAxis, const G4double pMin,const G4double pMax);
      // Restrict the volume to between specified min and max along the
      // given axis. Cartesian axes only, pMin<=pMax.

    G4double GetMaxXExtent() const;
      // Return maximum x extent.
    G4double GetMaxYExtent() const;
      // Return maximum y extent.
    G4double GetMaxZExtent() const;
      // Return maximum z extent.

    G4double GetMinXExtent() const;
      // Return minimum x extent.
    G4double GetMinYExtent() const;
      // Return minimum y extent.
    G4double GetMinZExtent() const;
      // Return minimum z extent.

    G4double GetMaxExtent(const EAxis pAxis) const;
      // Return maximum extent of volume along specified axis.
    G4double GetMinExtent(const EAxis pAxis) const;
      // Return minimum extent of volume along specified axis.

    G4bool IsXLimited() const;
      // Return true if the x axis is limited.
    G4bool IsYLimited() const;
      // Return true if the y axis is limited.
    G4bool IsZLimited() const;
      // Return true if the z axis is limited.

    G4bool IsLimited() const;
      // Return true if limited along any axis
    G4bool IsLimited(const EAxis pAxis) const;
      // Return true if the specified axis is restricted/limited.

    G4bool ClipToLimits(G4ThreeVector& pStart,G4ThreeVector& pEnd) const;
      // Clip the line segment pStart->pEnd to the volume described by the
      // current limits. Return true if the line remains after clipping,
      // else false, and leave the vectors in an undefined state.

    G4bool Inside(const G4ThreeVector& pVec) const;
      // Return true if the specified vector is inside/on boundaries of limits.

    G4int OutCode(const G4ThreeVector& pVec) const;
      // Calculate the `outcode' for the specified vector.
      // Intended for use during clipping against the limits
      // The bits are set given the following conditions:
      //   0      pVec.x()<fxAxisMin && IsXLimited()
      //   1      pVec.x()>fxAxisMax && IsXLimited()
      //   2      pVec.y()<fyAxisMin && IsYLimited()
      //   3      pVec.y()>fyAxisMax && IsYLimited()
      //   4      pVec.z()<fzAxisMin && IsZLimited()
      //   5      pVec.z()>fzAxisMax && IsZLimited()

  private:

    G4double fxAxisMin,fxAxisMax;
    G4double fyAxisMin,fyAxisMax;
    G4double fzAxisMin,fzAxisMax;
};

#include "G4VoxelLimits.icc"

G4std::ostream& operator << (G4std::ostream& os, const G4VoxelLimits& pLim);
  // Print the limits to the stream in the form:
  //  "{(xmin,xmax) (ymin,ymax) (zmin,zmax)}"
  // Replace (xmin,xmax) by (-,-)  when not limited.

#endif
