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
// G4Voxelizer
//
// Class description:
//
// Voxelizer for tessellated surfaces and solids positioning in 3D space,
// used in G4TessellatedSolid and G4MultiUnion.

// Author: Marek Gayer (CERN), 19.10.2012 - Created
// --------------------------------------------------------------------
#ifndef G4VOXELIZER_HH
#define G4VOXELIZER_HH

#include <vector>
#include <string>
#include <map>

#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4SurfBits.hh"
#include "G4Box.hh"
#include "G4VFacet.hh"
#include "G4VSolid.hh"

struct G4VoxelBox
{
  G4ThreeVector hlen; // half length of the box
  G4ThreeVector pos; // position of the box
};

struct G4VoxelInfo
{
  G4int count;
  G4int previous;
  G4int next;
};

/**
 * @brief G4Voxelizer is a tool for generating the optimisation structure
 * of tessellated surfaces and solids positioned in 3D space; it is used in
 * G4TessellatedSolid and G4MultiUnion.
 */

class G4Voxelizer
{
  public:

    /**
     * Constructor and default Destructor.
     */
    G4Voxelizer();
    ~G4Voxelizer() = default;

    /**
     * Builds the voxelisation structure for solids positioned in space.
     *  @param[in] solids The list of solids.
     *  @param[in] transforms The associated transformation in space.
     */
    void Voxelize(std::vector<G4VSolid*>& solids,
                  std::vector<G4Transform3D>& transforms);

    /**
     * Builds the voxelisation structure for facets forming a shape.
     *  @param[in] facets The list of facets.
     */
    void Voxelize(std::vector<G4VFacet*>& facets);

    /**
     * Displays the dX, dY, dZ, pX, pY and pZ for each node.
     */
    void DisplayVoxelLimits() const;

    /**
     * Prints the positions of the boundaries of the slices on the three axes.
     */
    void DisplayBoundaries();

    /**
     * Prints which solids are present in the slices previously elaborated.
     */
    void DisplayListNodes() const;

    /**
     * Displays the nodes located in a voxel characterised by its three indexes.
     */
    void GetCandidatesVoxel(std::vector<G4int>& voxels);

    /**
     * Methods returning in a vector container the nodes located in a voxel
     * characterised by its three indexes.
     *  @returns The total candidates number.
     */
    G4int GetCandidatesVoxelArray(const G4ThreeVector& point,
                                        std::vector<G4int>& list,
                                        G4SurfBits* crossed = nullptr) const;
    G4int GetCandidatesVoxelArray(const std::vector<G4int>& voxels,
                                  const G4SurfBits bitmasks[],
                                        std::vector<G4int>& list,
                                        G4SurfBits* crossed = nullptr) const;
    G4int GetCandidatesVoxelArray(const std::vector<G4int>& voxels,
                                        std::vector<G4int>& list,
                                        G4SurfBits* crossed = nullptr)const;

    /**
     * Returns the pointer to the array containing the characteristics
     * of each box.
     */
    inline const std::vector<G4VoxelBox>& GetBoxes() const;

    /**
     * Returns the boundary vector, given an 'index'.
     */
    inline const std::vector<G4double>& GetBoundary(G4int index) const;

    /**
     * Utility method for checking/updating current voxel given in input.
     */
    G4bool UpdateCurrentVoxel(const G4ThreeVector& point,
                              const G4ThreeVector& direction,
                                    std::vector<G4int>& curVoxel) const;

    /**
     * Updates current voxel based on provided 'point'.
     */
    inline void GetVoxel(std::vector<G4int>& curVoxel,
                         const G4ThreeVector& point) const;

    /**
     * Returns memory size of a slice.
     */
    inline G4int GetBitsPerSlice () const;

    /**
     * Returns true if 'point' is contained within boundaries.
     */
    G4bool Contains(const G4ThreeVector& point) const;

    /**
     * Returns the distance to next boundary, given 'point' and 'direction'
     * and updates current voxel 'curVoxel', using the index corresponding
     * to the closest voxel boundary on the ray.
     */
    G4double DistanceToNext(const G4ThreeVector& point,
                            const G4ThreeVector& direction,
                                  std::vector<G4int>& curVoxel) const;

    /**
     * Returns the distance to first bounding box, given 'point' and 'direction'.
     */
    G4double DistanceToFirst(const G4ThreeVector& point,
                             const G4ThreeVector& direction) const;

    /**
     * Returns the minimum distance of 'point' to the bounding box.
     */
    G4double DistanceToBoundingBox(const G4ThreeVector& point) const;

    /**
     * Utility for estimating the isotropic safety from a point 'p' outside
     * the current solid to any of its surfaces. The algorithm may be accurate
     * or should provide a fast underestimate, based on safety point 'f'.
     */
    G4double MinDistanceToBox (const G4ThreeVector& p,
                               const G4ThreeVector& f) const;

    /**
     * Accessors for voxels and points.
     */
    inline G4int GetVoxelsIndex(G4int x, G4int y, G4int z) const;
    inline G4int GetVoxelsIndex(const std::vector<G4int>& voxels) const;
    inline G4bool GetPointVoxel(const G4ThreeVector& p,
                                std::vector<G4int>& voxels) const;
    inline G4int GetPointIndex(const G4ThreeVector& p) const;

    /**
     * Returns the empty bits container.
     */
    inline const G4SurfBits& Empty() const;

    /**
     * Returns true if empty bit in container, given an 'index'.
     */
    inline G4bool IsEmpty(G4int index) const;

    /**
     * Setters/getter for the maximum number of voxels.
     */
    void SetMaxVoxels(G4int max);
    void SetMaxVoxels(const G4ThreeVector& reductionRatio);
    inline G4int GetMaxVoxels(G4ThreeVector& ratioOfReduction);

    /**
     * Logger returning the size of total allocated memory.
     */
    G4int AllocatedMemory();

    /**
     * Utility accessors/functions for voxels.
     */
    inline long long GetCountOfVoxels() const;
    inline long long CountVoxels(std::vector<G4double> boundaries[]) const;
    inline const std::vector<G4int>&
                 GetCandidates(std::vector<G4int>& curVoxel) const;
    inline G4int GetVoxelBoxesSize() const;
    inline const G4VoxelBox& GetVoxelBox(G4int i) const;
    inline const std::vector<G4int>& GetVoxelBoxCandidates(G4int i) const;
    inline G4int GetTotalCandidates() const;
    static void SetDefaultVoxelsCount(G4int count);
    static G4int GetDefaultVoxelsCount();

  private:

    /**
     * Binary search function for retrieving a value in a vector.
     */
    template <typename T> 
    inline G4int BinarySearch(const std::vector<T>& vec, T value) const;

    /**
     * Utilities.
     */
    G4String GetCandidatesAsString(const G4SurfBits& bits) const;
    void CreateSortedBoundary(std::vector<G4double>& boundaryRaw, G4int axis);
    void DisplayBoundaries(std::vector<G4double>& fBoundaries);
    void FindComponentsFastest(unsigned int mask,
                               std::vector<G4int>& list, G4int i) const;
    inline G4ThreeVector GetGlobalPoint(const G4Transform3D& trans,
                                        const G4ThreeVector& lpoint) const;
    void TransformLimits(G4ThreeVector& min, G4ThreeVector& max,
                         const G4Transform3D& transformation) const;

    /**
     * Build utilities.
     */
    void BuildEmpty ();
    void BuildBoundaries();
    void BuildReduceVoxels(std::vector<G4double> fBoundaries[],
                           G4ThreeVector reductionRatio);
    void BuildReduceVoxels2(std::vector<G4double> fBoundaries[],
                            G4ThreeVector reductionRatio);
    void BuildVoxelLimits(std::vector<G4VSolid*>& solids,
                          std::vector<G4Transform3D>& transforms);
    void BuildVoxelLimits(std::vector<G4VFacet*>& facets);
    void CreateMiniVoxels(std::vector<G4double> fBoundaries[],
                          G4SurfBits bitmasks[]);
    void BuildBitmasks(std::vector<G4double> fBoundaries[],
                       G4SurfBits bitmasks[], G4bool countsOnly = false);
    void BuildBoundingBox();
    void BuildBoundingBox(G4ThreeVector& amin, G4ThreeVector& amax,
                          G4double tolerance = 0.0);
    void SetReductionRatio(G4int maxVoxels, G4ThreeVector& reductionRatio);


  private:

   /**
    * @brief G4VoxelComparator is utility class used for comparing voxels.
    */

    class G4VoxelComparator
    {
      public:

      std::vector<G4VoxelInfo>& fVoxels;

      G4VoxelComparator(std::vector<G4VoxelInfo>& voxels) : fVoxels(voxels) {}

      G4bool operator()(const G4int& l, const G4int& r) const
      {
        G4VoxelInfo &lv = fVoxels[l], &rv = fVoxels[r];
        G4int left = lv.count +  fVoxels[lv.next].count;
        G4int right = rv.count + fVoxels[rv.next].count;
        return (left == right) ? l < r : left < right;
      }
    };

    static G4int fDefaultVoxelsCount;

    std::vector<G4VoxelBox> fVoxelBoxes;
    std::vector<std::vector<G4int> > fVoxelBoxesCandidates;
    mutable std::map<G4int, std::vector<G4int> > fCandidates;

    const std::vector<G4int> fNoCandidates;

    long long fCountOfVoxels;

    G4int fNPerSlice;

    std::vector<G4VoxelBox> fBoxes;
      // Array of box limits on the 3 cartesian axis

    std::vector<G4double> fBoundaries[3];
      // Sorted and if need skimmed fBoundaries along X,Y,Z axis

    std::vector<G4int> fCandidatesCounts[3]; 

    G4int fTotalCandidates;

    G4SurfBits fBitmasks[3];

    G4ThreeVector fBoundingBoxCenter;

    G4Box fBoundingBox;

    G4ThreeVector fBoundingBoxSize;

    G4ThreeVector fReductionRatio;

    G4int fMaxVoxels;

    G4double fTolerance;

    G4SurfBits fEmpty;
};

#include "G4Voxelizer.icc"

#endif
