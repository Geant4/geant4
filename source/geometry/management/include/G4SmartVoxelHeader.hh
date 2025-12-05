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
// G4SmartVoxelHeader
//
// Class description:
//
// Represents a set of voxels, created by a single axis of virtual division.
// Contains the individual voxels, which are potentially further divided
// along different axes.
//
// Member data:
//
// EAxis faxis
//   - The (cartesian) slicing/division axis
// G4double fmaxExtent
// G4double fminExtent
//   - Minimum and maximum coordiantes along the axis
// std::vector<G4SmartVoxelProxy*> fslices
//   - The slices along the axis
//
// G4int fminEquivalent
// G4int fmaxEquivalent
//   - Minimum and maximum equivalent slice nos.
//     [Applies to the level of the header, not its nodes]

// Author: Paul Kent (CERN), 13.07.1995 - Initial version
//         Gabriele Cosmo (CERN), 18.04.2001 - Migrated to STL vector
// --------------------------------------------------------------------
#ifndef G4SMARTVOXELHEADER_HH
#define G4SMARTVOXELHEADER_HH

#include <vector>

#include "G4Types.hh"
#include "geomdefs.hh"

#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"

// Forward declarations
class G4LogicalVolume;
class G4VoxelLimits;
class G4VPhysicalVolume;

// Aliases
using G4ProxyVector = std::vector<G4SmartVoxelProxy*>;
using G4NodeVector = std::vector<G4SmartVoxelNode*>;
using G4VolumeNosVector = std::vector<G4int>;
using G4VolumeExtentVector = std::vector<G4double>;

/**
 * @brief G4SmartVoxelHeader represents a set of voxels, created by a single
 * axis of virtual division. It contains the individual voxels, which are
 * potentially further divided along different axes.
 */

class G4SmartVoxelHeader
{
  public:

    /**
     * Constructor for the topmost header, to begin voxel construction at a
     * given logical volume; 'pSlice' is used to set max and min equivalent
     * slice numbers for the header - they apply to the level of the header,
     * not its nodes.
     *  @param[in] pVolume Pointer to the logical volume to voxelise.
     *  @param[in] pSlice Max & min equivalent slice numbers for the header.
     */
    G4SmartVoxelHeader(G4LogicalVolume* pVolume, G4int pSlice = 0);

    /**
     * Constructor for building and refine voxels between specified limits,
     * considering only the physical volumes numbered 'pCandidates'; 'pSlice'
     * is used to set max and min equivalent slice numbers for the header;
     * they apply to the level of the header, not its nodes.
     *  @param[in] pVolume Pointer to the logical volume to voxelise.
     *  @param[in] pLimits Refinement limits for building the voxels.
     *  @param[in] pCandidates Candidate volumes to be considered.
     *  @param[in] pSlice Max & min equivalent slice numbers for the header.
     */
    G4SmartVoxelHeader(G4LogicalVolume* pVolume,
		       const G4VoxelLimits& pLimits,
		       const G4VolumeNosVector* pCandidates,
		       G4int pSlice = 0);
      // 

    /**
     * Destructor. Deletes all referenced nodes [but *not* the referenced
     * physical volumes].
     */
    ~G4SmartVoxelHeader();
    
    /**
     * Equality operator.
     */
    G4bool operator == (const G4SmartVoxelHeader& pHead) const;

    /**
     * Streaming operator.
     */
    friend std::ostream&
    operator << (std::ostream&s, const G4SmartVoxelHeader& h);

    /**
     * Access functions for min/max equivalent slices (nodes & headers).
     */
    G4int GetMaxEquivalentSliceNo() const;
    void SetMaxEquivalentSliceNo(G4int pMax);
    G4int GetMinEquivalentSliceNo() const;
    void SetMinEquivalentSliceNo(G4int pMin);

    /**
     * Returns the current division axis.
     */
    EAxis GetAxis() const;

    /**
     * Returns the suggested division axis for parameterised volume.
     */
    EAxis GetParamAxis() const;

    /**
     * Returns the maximum coordinate limit along the current axis.
     */
    G4double GetMaxExtent() const;

    /**
     * Returns the minimum coordinate limit along the current axis.
     */
    G4double GetMinExtent() const;
    
    /**
     * Returns the number of slices along the current axis.
     */
    std::size_t GetNoSlices() const;
    
    /**
     * Returns the pointer to the proxy for the n-th slice
     * (numbering from 0, no bounds checking is performed).
     */
    G4SmartVoxelProxy* GetSlice(std::size_t n) const;

    /**
     * Returns true if all slices are equal (after collection).
     */
    G4bool AllSlicesEqual() const;

  private:  // 'Worker' / operation functions

    /**
     * Builds and refine voxels for daughters of a specified volume 'pVolume'
     * which DOES NOT contain a REPLICATED daughter.
     */
    void BuildVoxels(G4LogicalVolume* pVolume);

    /**
     * Builds voxels for a specified volume 'pVolume' containing a single
     * replicated volume.
     */
    void BuildReplicaVoxels(G4LogicalVolume* pVolume);

    /**
     * Constructs nodes in simple consuming case.
     */
    void BuildConsumedNodes(G4int nReplicas);

    /**
     * Builds and refines voxels between specified limits 'pLimits',
     * considering only the physical volumes 'pCandidates'. Main entry point
     * for "construction". Hardwired to stop at third level of refinement,
     * using the XYZ Cartesian axes in any order.
     */
    void BuildVoxelsWithinLimits(G4LogicalVolume* pVolume,
                                 G4VoxelLimits pLimits,
		                 const G4VolumeNosVector* pCandidates);

    /**
     * Calculates and stores the minimum and maximum equivalent neighbour
     * values for all slices.
     */
    void BuildEquivalentSliceNos();

    /**
     * Collects the common nodes, deleting all but one to save memory
     * and adjusting stored slice pointers appropriately.
     */
    void CollectEquivalentNodes();

    /**
     * Collects the common headers, deleting all but one to save memory
     * and adjusting stored slice pointers appropriately.
     */
    void CollectEquivalentHeaders();

    /**
     * Builds the nodes corresponding to the specified axis 'pAxis', within
     * the specified limits 'pLimits', considering the daughters numbered
     * 'pCandidates' of the logical volume 'pVolume'.
     */
    G4ProxyVector* BuildNodes(G4LogicalVolume* pVolume,
	                      G4VoxelLimits pLimits,
	                      const G4VolumeNosVector* pCandidates,
	                      EAxis pAxis);

    /**
     * Calculates a "quality value" for the specified vector of voxels.
     * The value returned should be greater than zero and such that the
     * smaller the number the higher the quality of the slice; 'pSlice'
     * must consist of smart voxel node proxies only.
     */
    G4double CalculateQuality(G4ProxyVector* pSlice);

    /**
     * Examined each contained node, refines (creates a replacement additional
     * dimension of voxels) when there is more than one voxel in the slice.
     */
    void RefineNodes(G4LogicalVolume* pVolume, G4VoxelLimits pLimits);

  private:

    /** Min and max equivalent slice nos for previous level. */
    G4int fminEquivalent;
    G4int fmaxEquivalent;

    /** Axis for slices. */
    EAxis faxis, fparamAxis;

    /** Max and min coordinate along faxis. */
    G4double fmaxExtent;
    G4double fminExtent;

    /** Slices along axis. */
    G4ProxyVector fslices;
};

#include "G4SmartVoxelHeader.icc"

#endif
