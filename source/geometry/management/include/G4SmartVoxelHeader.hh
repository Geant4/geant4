// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SmartVoxelHeader.hh,v 1.4 2000-04-20 16:49:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SmartVoxelHeader
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
// G4RWTPtrOrderedVector<G4SmartVoxelProxy> fslices
//   - The slices along the axis
//
// G4int fminEquivalent
// G4int fmaxEquivalent
//   - Minimum and maximum equivalent slice nos.
//     [Applies to the level of the header, not its nodes]

// History:
// 13.07.95 P.Kent Initial version

#ifndef G4SMARTVOXELHEADER_HH
#define G4SMARTVOXELHEADER_HH

#include "globals.hh"
#include "geomdefs.hh"
#include "voxeldefs.hh"

#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"

#include "G4ios.hh"

#include "g4rw/tvvector.h"
#include "g4rw/tpordvec.h"
#include "g4rw/tvordvec.h"

// Forward declarations
class G4LogicalVolume;
class G4VoxelLimits;
class G4VPhysicalVolume;

// Typedefs
typedef G4RWTPtrOrderedVector<G4SmartVoxelProxy> G4ProxyVector;
typedef G4RWTPtrOrderedVector<G4SmartVoxelNode> G4NodeVector;
typedef G4RWTValOrderedVector<G4int> G4VolumeNosVector;
typedef G4RWTValVector<G4double> G4VolumeExtentVector;

class G4SmartVoxelHeader
{
  public:  // with description

    G4SmartVoxelHeader(G4LogicalVolume* pVolume,const G4int pSlice=0);
      // Constructor for topmost header, to begin voxel construction at a
      // given logical volume. pSlice is used to set max and min equivalent
      // slice nos for the header - they apply to the level of the header,
      // not its nodes.

    ~G4SmartVoxelHeader();
      // Delete all referenced nodes [but *not* referenced physical volumes].
    
    G4int GetMaxEquivalentSliceNo() const;
    void SetMaxEquivalentSliceNo(const G4int pMax);
    G4int GetMinEquivalentSliceNo() const;
    void SetMinEquivalentSliceNo(const G4int pMin);
      // Access functions for min/max equivalent slices (nodes & headers).

    EAxis GetAxis() const;
      // Return the current division axis.
    
    G4double GetMaxExtent() const;
      // Return the maximum coordinate limit along the current axis.
    G4double GetMinExtent() const;
      // Return the minimum coordinate limit along the current axis.
    
    G4int GetNoSlices() const;
      // Return the no of slices along the current axis.
    
    G4SmartVoxelProxy* GetSlice(const G4int n) const;
      // Return ptr to the proxy for the nth slice (numbering from 0,
      // no bounds checking performed).

    G4bool AllSlicesEqual() const;
      // True if all slices equal (after collection).

  public:  // without description

    G4bool operator == (const G4SmartVoxelHeader& pHead) const;

    friend G4std::ostream&
    operator << (G4std::ostream&s, const G4SmartVoxelHeader& h);

  protected:

    G4SmartVoxelHeader(G4LogicalVolume* pVolume,
		       const G4VoxelLimits& pLimits,
		       const G4VolumeNosVector* pCandidates,
		       const G4int pSlice=0);
      // Build and refine voxels between specified limits, considering only
      // the physical volumes numbered `pCandidates'. pSlice is used to set max
      // and min equivalent slice nos for the header - they apply to the level
      // of the header, not its nodes.

    //  `Worker' / operation functions:

    void BuildVoxels(G4LogicalVolume* pVolume);
      // Build and refine voxels for daughters of specified volume which
      // DOES NOT contain a REPLICATED daughter.

    void BuildReplicaVoxels(G4LogicalVolume* pVolume);
      // Build voxels for specified volume containing a single
      // replicated volume.

    void BuildConsumedNodes(const G4int nReplicas);
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

    EAxis faxis;
      // Axis for slices.

    G4double fmaxExtent;
    G4double fminExtent;
      // Max and min coordinate along faxis.

    G4ProxyVector fslices;
      // Slices along axis.
};

#include "G4SmartVoxelHeader.icc"

#endif
