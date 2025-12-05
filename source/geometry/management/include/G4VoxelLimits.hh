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
// G4VoxelLimits
//
// Class description:
//
// Represents limitation/restrictions of space, where restrictions
// are only made perpendicular to the cartesian axes.
//
// Member data:
//
// G4double fxAxisMin,fxAxisMax
// G4double fyAxisMin,fyAxisMax
// G4double fzAxisMin,fzAxisMax
// - The min and max values along each axis. +-kInfinity if not restricted.

// Author: Paul Kent (CERN), 13.07.1995 - Initial version.
// --------------------------------------------------------------------
#ifndef G4VOXELLIMITS_HH
#define G4VOXELLIMITS_HH

#include "G4Types.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"

#include <assert.h>

/**
 * @brief G4VoxelLimits represents limitation/restrictions of space, where
 * restrictions are only made perpendicular to the Cartesian axes.
 */

class G4VoxelLimits
{
  public:
  
    /**
     * Default Constructor & Destructor.
     * Constructor initialises to be unlimited. Volume unrestricted.
     */
    G4VoxelLimits() = default;
    ~G4VoxelLimits() = default;

    /**
     * Restricts the volume to between specified min and max along the
     * given axis. Cartesian axes only, pMin<=pMax.
     */
    void AddLimit(const EAxis pAxis, const G4double pMin, const G4double pMax);

    /**
     * Accessors for max extent
     */
    G4double GetMaxXExtent() const;
    G4double GetMaxYExtent() const;
    G4double GetMaxZExtent() const;

    /**
     * Accessors for min extent
     */
    G4double GetMinXExtent() const;
    G4double GetMinYExtent() const;
    G4double GetMinZExtent() const;

    /**
     * Accessors for the extent of the volume along the specified axis.
     */
    G4double GetMaxExtent(const EAxis pAxis) const;
    G4double GetMinExtent(const EAxis pAxis) const;

    /**
     * Return true if the X/Y/Z axis is limited.
     */
    G4bool IsXLimited() const;
    G4bool IsYLimited() const;
    G4bool IsZLimited() const;

    /**
     * Return true if limited along any axis.
     */
    G4bool IsLimited() const;

    /**
     * Return true if the specified axis is restricted/limited.
     */
    G4bool IsLimited(const EAxis pAxis) const;

    /**
     * Clips the line segment pStart->pEnd to the volume described by the
     * current limits. Returns true if the line remains after clipping,
     * else false, and leaves the vectors in an undefined state.
     */
    G4bool ClipToLimits(G4ThreeVector& pStart, G4ThreeVector& pEnd) const;

    /**
     * Returns true if the specified vector is inside/on boundaries of limits.
     */
    G4bool Inside(const G4ThreeVector& pVec) const;

    /**
     * Calculates the 'outcode' for the specified vector.
     * Intended for use during clipping against the limits.
     * The bits are set given the following conditions:
     *   0      pVec.x()<fxAxisMin && IsXLimited()
     *   1      pVec.x()>fxAxisMax && IsXLimited()
     *   2      pVec.y()<fyAxisMin && IsYLimited()
     *   3      pVec.y()>fyAxisMax && IsYLimited()
     *   4      pVec.z()<fzAxisMin && IsZLimited()
     *   5      pVec.z()>fzAxisMax && IsZLimited()
     */
    G4int OutCode(const G4ThreeVector& pVec) const;

  private:

    G4double fxAxisMin = -kInfinity, fxAxisMax = kInfinity;
    G4double fyAxisMin = -kInfinity, fyAxisMax = kInfinity;
    G4double fzAxisMin = -kInfinity, fzAxisMax = kInfinity;
};

/**
 * Prints the limits to the stream in the form:
 * "{(xmin,xmax) (ymin,ymax) (zmin,zmax)}"
 * Replaces (xmin,xmax) by (-,-)  when not limited.
 */
std::ostream& operator << (std::ostream& os, const G4VoxelLimits& pLim);

#include "G4VoxelLimits.icc"

#endif
