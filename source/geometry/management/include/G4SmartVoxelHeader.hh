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

// 18.04.01, G.Cosmo - Migrated to STL vector
// 13.07.95, P.Kent - Initial version
// --------------------------------------------------------------------
#ifndef G4SMARTVOXELHEADER_HH
#define G4SMARTVOXELHEADER_HH 1

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

class G4SmartVoxelHeader
{
  public:

    G4SmartVoxelHeader(G4LogicalVolume* pVolume, G4int pSlice = 0);
      // Constructor for topmost header, to begin voxel construction at a
      // given logical volume. pSlice is used to set max and min equivalent
      // slice nos for the header - they apply to the level of the header,
      // not its nodes.

    ~G4SmartVoxelHeader();
      // Delete all referenced nodes [but *not* referenced physical volumes].
    
    G4int GetMaxEquivalentSliceNo() const;
    void SetMaxEquivalentSliceNo(G4int pMax);
    G4int GetMinEquivalentSliceNo() const;
    void SetMinEquivalentSliceNo(G4int pMin);
      // Access functions for min/max equivalent slices (nodes & headers).

    EAxis GetAxis() const;
      // Return the current division axis.
    EAxis GetParamAxis() const;
      // Return suggested division axis for parameterised volume.

    G4double GetMaxExtent() const;
      // Return the maximum coordinate limit along the current axis.
    G4double GetMinExtent() const;
      // Return the minimum coordinate limit along the current axis.
    
    std::size_t GetNoSlices() const;
      // Return the no of slices along the current axis.
    
    G4SmartVoxelProxy* GetSlice(std::size_t n) const;
      // Return ptr to the proxy for the nth slice (numbering from 0,
      // no bounds checking performed).

    G4bool AllSlicesEqual() const;
      // True if all slices equal (after collection).

  public:

    G4bool operator == (const G4SmartVoxelHeader& pHead) const;

    friend std::ostream&
    operator << (std::ostream&s, const G4SmartVoxelHeader& h);

    G4SmartVoxelHeader(G4LogicalVolume* pVolume,
		       const G4VoxelLimits& pLimits,
		       const G4VolumeNosVector* pCandidates,
		       G4int pSlice = 0);
      // Build and refine voxels between specified limits, considering only
      // the physical volumes numbered `pCandidates'. pSlice is used to set max
      // and min equivalent slice nos for the header - they apply to the level
      // of the header, not its nodes.

  protected:

    //  `Worker' / operation functions:

    void BuildVoxels(G4LogicalVolume* pVolume);
      // Build and refine voxels for daughters of specified volume which
      // DOES NOT contain a REPLICATED daughter.

    void BuildReplicaVoxels(G4LogicalVolume* pVolume);
      // Build voxels for specified volume containing a single
      // replicated volume.

    void BuildConsumedNodes(G4int nReplicas);
      // Construct nodes in simple consuming case.

    void BuildVoxelsWithinLimits(G4LogicalVolume* pVolume,
                                 G4VoxelLimits pLimits,
		                 const G4VolumeNosVector* pCandidates);
      // Build and refine voxels between specified limits, considering only
      // the physical volumes `pCandidates'. Main entry point for "construction".
      // Hardwired to stop at third level of refinement, using the xyz cartesian
      // axes in any order.

    void BuildEquivalentSliceNos();
      // Calculate and Store the minimum and maximum equivalent neighbour
      // values for all slices.

    void CollectEquivalentNodes();
      // Collect common nodes, deleting all but one to save memory,
      // and adjusting stored slice ptrs appropriately.

    void CollectEquivalentHeaders();
      // Collect common headers, deleting all but one to save memory,
      // and adjusting stored slice ptrs appropriately.


    G4ProxyVector* BuildNodes(G4LogicalVolume* pVolume,
	                      G4VoxelLimits pLimits,
	                      const G4VolumeNosVector* pCandidates,
	                      EAxis pAxis);
      // Build the nodes corresponding to the specified axis, within
      // the specified limits, considering the daughters numbered pCandidates
      // of the logical volume.

    G4double CalculateQuality(G4ProxyVector *pSlice);
      // Calculate a "quality value" for the specified vector of voxels
      // The value returned should be >0 and such that the smaller the
      // number the higher the quality of the slice.
      // pSlice must consist of smartvoxelnodeproxies only.

    void RefineNodes(G4LogicalVolume* pVolume,G4VoxelLimits pLimits);
      // Examined each contained node, refine (create a replacement additional
      // dimension of voxels) when there is more than one voxel in the slice.

    G4int fminEquivalent;
    G4int fmaxEquivalent;
      // Min and max equivalent slice nos for previous level.

    EAxis faxis, fparamAxis;
      // Axis for slices.

    G4double fmaxExtent;
    G4double fminExtent;
      // Max and min coordinate along faxis.

    G4ProxyVector fslices;
      // Slices along axis.
};

#include "G4SmartVoxelHeader.icc"

#endif
