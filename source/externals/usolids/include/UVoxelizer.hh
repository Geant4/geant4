//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UVoxelizer
//
// Class description:
//
//    Voxelizer used for UPolycone, UPolyhedra, UTessellatedSolid
//    and UMultiUnion.
//
// 19.10.12 Marek Gayer
//          Created from original implementation in ROOT
// --------------------------------------------------------------------

#ifndef UVoxelizer_HH
#define UVoxelizer_HH

#include <vector>
#include <string>
#include <map>

#include "UBits.hh"
#include "UBox.hh"
#include "VUFacet.hh"
#include "VUSolid.hh"
#include "UUtils.hh"
#include "UTransform3D.hh"

struct UVoxelBox
{
  UVector3 hlen; // half length of the box
  UVector3 pos; // position of the box
};

struct UVoxelInfo
{
  int count;
  int previous;
  int next;
};

class UVoxelizer
{
//  friend class UVoxelCandidatesIterator;

  public:

    // Binary search
    template <typename T>
    static inline int BinarySearch(const std::vector<T>& vec, T value)
    {
      // Binary search in an array of doubles. If match is found, function returns
      // position of element.  If no match found, function gives nearest
      // element smaller than value.
      typename std::vector<T>::const_iterator begin = vec.begin(), end = vec.end();
      int res = std::upper_bound(begin, end, value) - begin - 1;
      return res;
    }

#ifdef USOLIDSONLY
    void Voxelize(std::vector<VUSolid*>& solids, std::vector<UTransform3D>& transforms);
#endif // USOLIDSONLY

    void Voxelize(std::vector<VUFacet*>& facets);

    void DisplayVoxelLimits();
    void DisplayBoundaries();
    void DisplayListNodes();

    UVoxelizer();
    ~UVoxelizer();

    // Method displaying the nodes located in a voxel characterized by its three indexes:
    void               GetCandidatesVoxel(std::vector<int>& voxels);
    // Method returning in a vector container the nodes located in a voxel characterized by its three indexes:
    int GetCandidatesVoxelArray(const UVector3& point, std::vector<int>& list, UBits* crossed = NULL) const;

    int GetCandidatesVoxelArray(const std::vector<int>& voxels, const UBits bitmasks[], std::vector<int>& list, UBits* crossed = NULL) const;

    int GetCandidatesVoxelArray(const std::vector<int>& voxels, std::vector<int>& list, UBits* crossed = NULL)const;

    // Method returning the pointer to the array containing the characteristics of each box:
    inline const std::vector<UVoxelBox>& GetBoxes() const
    {
      return fBoxes;
    }
    inline const std::vector<double>& GetBoundary(int index) const
    {
      return fBoundaries[index];
    }

    bool UpdateCurrentVoxel(const UVector3& point, const UVector3& direction, std::vector<int>& curVoxel) const;

    inline void GetVoxel(std::vector<int>& curVoxel, const UVector3& point) const
    {
      for (int i = 0; i <= 2; ++i)
      {
        const std::vector<double>& boundary = GetBoundary(i);
        int n = BinarySearch(boundary, point[i]);
        if (n == -1) n = 0;
        else if (n == (int) boundary.size() - 1) n--;
        curVoxel[i] = n;
      }
    }

    inline int GetBitsPerSlice() const
    {
      return fNPerSlice * 8 * sizeof(unsigned int);
    }

    bool Contains(const UVector3& point) const;

    double DistanceToNext(const UVector3& point, const UVector3& direction, std::vector<int>& curVoxel) const;

    double DistanceToFirst(const UVector3& point, const UVector3& direction) const;

    double SafetyToBoundingBox(const UVector3& point) const;

    inline int GetVoxelsIndex(int x, int y, int z) const
    {
      if (x < 0 || y < 0 || z < 0) return -1;
      int maxX = fBoundaries[0].size();
      int maxY = fBoundaries[1].size();
      int index = x + y * maxX + z * maxX * maxY;
      return index;
    }

    inline int GetVoxelsIndex(const std::vector<int>& voxels) const
    {
      return GetVoxelsIndex(voxels[0], voxels[1], voxels[2]);
    }

    inline bool GetPointVoxel(const UVector3& p, std::vector<int>& voxels) const
    {
      for (int i = 0; i <= 2; ++i)
        if (p[i] < *fBoundaries[i].begin() || p[i] > *fBoundaries[i].end()) return false;

      for (int i = 0; i <= 2; ++i)
        voxels[i] = BinarySearch(fBoundaries[i], p[i]);

      return true;
    }

    inline int GetPointIndex(const UVector3& p) const
    {
      int maxX = fBoundaries[0].size();
      int maxY = fBoundaries[1].size();
      int x = BinarySearch(fBoundaries[0], p[0]);
      int y = BinarySearch(fBoundaries[1], p[1]);
      int z = BinarySearch(fBoundaries[2], p[2]);
      int index = x + y * maxX + z * maxX * maxY;
      return index;
    }

    inline const UBits& Empty() const
    {
      return fEmpty;
    }

    inline bool IsEmpty(int index) const
    {
      return fEmpty[index];
    }

    void SetMaxVoxels(int max);

    void SetMaxVoxels(const UVector3& reductionRatio);

    inline int GetMaxVoxels(UVector3& ratioOfReduction)
    {
      ratioOfReduction = fReductionRatio;
      return fMaxVoxels;
    }

    int AllocatedMemory();

    inline long long GetCountOfVoxels() const
    {
      return fCountOfVoxels;
    }

    inline long long CountVoxels(std::vector<double> boundaries[]) const
    {
      long long sx = boundaries[0].size() - 1;
      long long sy = boundaries[1].size() - 1;
      long long sz = boundaries[2].size() - 1;
      return  sx * sy * sz;
    }

    inline const std::vector<int>& GetCandidates(std::vector<int>& curVoxel) const
    {
      int voxelsIndex = GetVoxelsIndex(curVoxel);
      if (voxelsIndex >= 0 && !fEmpty[voxelsIndex])
      {
        return fCandidates[voxelsIndex];
      }
      return fNoCandidates;
    }

    inline int GetVoxelBoxesSize() const
    {
      return fVoxelBoxes.size();
    }

    inline const UVoxelBox& GetVoxelBox(int i) const
    {
      return fVoxelBoxes[i];
    }

    inline const std::vector<int>& GetVoxelBoxCandidates(int i) const
    {
      return fVoxelBoxesCandidates[i];
    }

    inline int GetTotalCandidates() const
    {
      return fTotalCandidates;
    }

    static double MinDistanceToBox(const UVector3& aPoint, const UVector3& f);

    static void SetDefaultVoxelsCount(int count);

    static int GetDefaultVoxelsCount();

    void BuildBoundingBox();

    void BuildBoundingBox(UVector3& amin, UVector3& amax, double tolerance = 0);

    static void FindComponentsFastest(unsigned int mask,
                                      std::vector<int> &list, int i);
  private:

    static int fDefaultVoxelsCount;

    std::vector<UVoxelBox> fVoxelBoxes;

    std::vector<std::vector<int> > fVoxelBoxesCandidates;

    mutable std::map<int, std::vector<int> > fCandidates;

    const std::vector<int> fNoCandidates;

    long long fCountOfVoxels;

    void BuildEmpty();

    std::string GetCandidatesAsString(const UBits& bits);

    void CreateSortedBoundary(std::vector<double>& boundaryRaw, int axis);

    void BuildBoundaries();

    void BuildReduceVoxels(std::vector<double> fBoundaries[], UVector3 reductionRatio);

    void BuildReduceVoxels2(std::vector<double> fBoundaries[], UVector3 reductionRatio);

#ifdef USOLIDSONLY
    void BuildVoxelLimits(std::vector<VUSolid*>& solids, std::vector<UTransform3D>& transforms);
#endif // USOLIDSONLY

    void BuildVoxelLimits(std::vector<VUFacet*>& facets);

    void DisplayBoundaries(std::vector<double>& fBoundaries);

    void BuildBitmasks(std::vector<double> fBoundaries[], UBits bitmasks[]);

    void SetReductionRatio(int maxVoxels, UVector3& reductionRatio);

    void CreateMiniVoxels(std::vector<double> fBoundaries[], UBits bitmasks[]);

    int fNPerSlice;

    std::vector<UVoxelBox> fBoxes; // Array of box limits on the 3 cartesian axis

    std::vector<double> fBoundaries[3]; // Sorted and if need skimmed fBoundaries along X,Y,Z axis

    std::vector<int> fCandidatesCounts[3];

    int fTotalCandidates;

    UBits fBitmasks[3];

    UVector3 fBoundingBoxCenter;
    UBox fBoundingBox;
    UVector3 fBoundingBoxSize;

    UVector3 fReductionRatio;

    int fMaxVoxels;

    double fTolerance;

    UBits fEmpty;
};

#endif
