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
// $Id:$
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Voxelizer
//
// Class description:
//
// Voxelizer for tessellated surfaces and solids positioning in 3D space,
// used in G4TessellatedSolid and G4MultiUnion.

// History:
// 19.10.12 Marek Gayer, created
// --------------------------------------------------------------------

#ifndef G4Voxelizer_HH
#define G4Voxelizer_HH

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

class G4Voxelizer
{
  public:

    template <typename T> 
    static inline G4int BinarySearch(const std::vector<T> &vec, T value);

    void Voxelize(std::vector<G4VSolid *> &solids,
                  std::vector<G4Transform3D>& transforms);
    void Voxelize(std::vector<G4VFacet *> &facets);

    void DisplayVoxelLimits() const;
    void DisplayBoundaries();
    void DisplayListNodes() const;

    G4Voxelizer();
   ~G4Voxelizer();

    void GetCandidatesVoxel(std::vector<G4int> &voxels);
      // Method displaying the nodes located in a voxel characterized
      // by its three indexes.

    G4int GetCandidatesVoxelArray(const G4ThreeVector &point,
                                        std::vector<G4int> &list,
                                        G4SurfBits *crossed=0) const;
      // Method returning in a vector container the nodes located in a voxel
      // characterized by its three indexes.
    G4int GetCandidatesVoxelArray(const std::vector<G4int> &voxels,
                                  const G4SurfBits bitmasks[],
                                        std::vector<G4int> &list,
                                        G4SurfBits *crossed=0) const;
    G4int GetCandidatesVoxelArray(const std::vector<G4int> &voxels,
                                        std::vector<G4int> &list,
                                        G4SurfBits *crossed=0)const;

    inline const std::vector<G4VoxelBox> &GetBoxes() const;
      // Method returning the pointer to the array containing the
      // characteristics of each box.

    inline const std::vector<G4double> &GetBoundary(G4int index) const;

    G4bool UpdateCurrentVoxel(const G4ThreeVector &point,
                              const G4ThreeVector &direction,
                                    std::vector<G4int> &curVoxel) const;

    inline void GetVoxel(std::vector<G4int> &curVoxel,
                         const G4ThreeVector &point) const;

    inline G4int GetBitsPerSlice () const;

    G4bool Contains(const G4ThreeVector &point) const;

    G4double DistanceToNext(const G4ThreeVector &point,
                            const G4ThreeVector &direction,
                                  std::vector<G4int> &curVoxel) const;

    G4double DistanceToFirst(const G4ThreeVector &point,
                             const G4ThreeVector &direction) const;

    G4double DistanceToBoundingBox(const G4ThreeVector &point) const;

    inline G4int GetVoxelsIndex(G4int x, G4int y, G4int z) const;
    inline G4int GetVoxelsIndex(const std::vector<G4int> &voxels) const;
    inline G4bool GetPointVoxel(const G4ThreeVector& p,
                                std::vector<G4int>& voxels) const;
    inline G4int GetPointIndex(const G4ThreeVector &p) const;

    inline const G4SurfBits &Empty() const;
    inline G4bool IsEmpty(G4int index) const;

    void SetMaxVoxels(G4int max);
    void SetMaxVoxels(const G4ThreeVector &reductionRatio);

    inline G4int GetMaxVoxels(G4ThreeVector &ratioOfReduction);

    G4int AllocatedMemory();

    inline long long GetCountOfVoxels() const;

    inline long long CountVoxels(std::vector<G4double> boundaries[]) const;

    inline const std::vector<G4int> &
                 GetCandidates(std::vector<G4int> &curVoxel) const;

    inline G4int GetVoxelBoxesSize() const;

    inline const G4VoxelBox &GetVoxelBox(G4int i) const;

    inline const std::vector<G4int> &GetVoxelBoxCandidates(G4int i) const;

    inline G4int GetTotalCandidates() const;

    static G4double MinDistanceToBox (const G4ThreeVector &aPoint,
                                      const G4ThreeVector &f);

    static void SetDefaultVoxelsCount(G4int count);

    static G4int GetDefaultVoxelsCount();

  private:

    class G4VoxelComparator
    {
      public:

      std::vector<G4VoxelInfo> &fVoxels;

      G4VoxelComparator(std::vector<G4VoxelInfo> &voxels) : fVoxels(voxels) {}

      G4bool operator()(const G4int& l, const G4int& r) const
      {
        G4VoxelInfo &lv = fVoxels[l], &rv = fVoxels[r];
        G4int left = lv.count +  fVoxels[lv.next].count;
        G4int right = rv.count + fVoxels[rv.next].count;
        return (left == right) ? l < r : left < right;
      }
    };

  private:

    void BuildEmpty ();

    G4String GetCandidatesAsString(const G4SurfBits &bits) const;

    void CreateSortedBoundary(std::vector<G4double> &boundaryRaw, G4int axis);

    void BuildBoundaries();

    void BuildReduceVoxels(std::vector<G4double> fBoundaries[],
                           G4ThreeVector reductionRatio);
    void BuildReduceVoxels2(std::vector<G4double> fBoundaries[],
                            G4ThreeVector reductionRatio);

    void BuildVoxelLimits(std::vector<G4VSolid*>& solids,
                          std::vector<G4Transform3D>& transforms);
    void BuildVoxelLimits(std::vector<G4VFacet *> &facets);

    void DisplayBoundaries(std::vector<G4double> &fBoundaries);

    void BuildBitmasks(std::vector<G4double> fBoundaries[],
                       G4SurfBits bitmasks[], G4bool countsOnly=false);

    void BuildBoundingBox();
    void BuildBoundingBox(G4ThreeVector& amin, G4ThreeVector& amax,
                          G4double tolerance = 0);

    void SetReductionRatio(G4int maxVoxels, G4ThreeVector &reductionRatio);

    void CreateMiniVoxels(std::vector<G4double> fBoundaries[],
                          G4SurfBits bitmasks[]);
    static void FindComponentsFastest(unsigned int mask,
                                      std::vector<G4int> &list, G4int i);

    inline G4ThreeVector GetGlobalPoint(const G4Transform3D& trans,
                                        const G4ThreeVector& lpoint) const;
    void TransformLimits(G4ThreeVector& min, G4ThreeVector& max,
                         const G4Transform3D& transformation) const;

  private:

    static G4ThreadLocal G4int fDefaultVoxelsCount;

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
