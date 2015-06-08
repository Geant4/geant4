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
// 19.10.12 Marek Gayer
//          Created from original implementation in ROOT
// --------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <set>

#include "VUSolid.hh"
#include "UUtils.hh"
#include "UOrb.hh"
#include "UVoxelizer.hh"

using namespace std;

//______________________________________________________________________________
UVoxelizer::UVoxelizer() : fBoundingBox("TessBBox", 1, 1, 1)
{
  fCountOfVoxels = fNPerSlice = fTotalCandidates = 0;

  fTolerance = VUSolid::Tolerance();

  SetMaxVoxels(fDefaultVoxelsCount);
}

//______________________________________________________________________________
UVoxelizer::~UVoxelizer()
{
}


void UVoxelizer::BuildEmpty()
{
  // by reserving the size of candidates, we would avoid reallocation of
  // the vector which could cause fragmentation
  vector<int> xyz(3), candidates(fTotalCandidates);
  const vector<int> empty(0);
  int max[3];

  for (int i = 0; i <= 2; ++i) max[i] = fBoundaries[i].size();
  unsigned int size = max[0] * max[1] * max[2];

  fEmpty.Clear();
  fEmpty.ResetBitNumber(size - 1);
  fEmpty.ResetAllBits(true);

  for (xyz[2] = 0; xyz[2] < max[2]; ++xyz[2])
  {
    for (xyz[1] = 0; xyz[1] < max[1]; ++xyz[1])
    {
      for (xyz[0] = 0; xyz[0] < max[0]; ++xyz[0])
      {
        if (GetCandidatesVoxelArray(xyz, candidates))
        {
          int index = GetVoxelsIndex(xyz);
          fEmpty.SetBitNumber(index, false);
          // rather than assigning directly with:
          // "fCandidates[index] = candidates;", in an effort to ensure that
          // capacity would be just exact, we rather use following 3 lines:
          vector<int>& c = (fCandidates[index] = empty);
          c.reserve(candidates.size());
          c.assign(candidates.begin(), candidates.end());
        }
      }
    }
  }
#ifdef USPECSDEBUG
  cout << "Non-empty voxels count: " << fCandidates.size() << endl;
#endif
}

#ifdef USOLIDSONLY
//______________________________________________________________________________
//void UVoxelizer::BuildVoxelLimits(vector<VUSolid*>& solids, vector<UTransform3D*>& transforms)
void UVoxelizer::BuildVoxelLimits(vector<VUSolid*>& solids, vector<UTransform3D>& transforms)
{
  // "BuildVoxelLimits"'s aim is to store the coordinates of the origin as well as
  // the half lengths related to the bounding box of each node.
  // These quantities are stored in the array "fBoxes" (6 different values per
  // node.
  if (int numNodes = solids.size()) // Number of nodes in "multiUnion"
  {
    fBoxes.resize(numNodes); // Array which will store the half lengths
    fNPerSlice = 1 + (fBoxes.size() - 1) / (8 * sizeof(unsigned int));

    // related to a particular node, but also
    // the coordinates of its origin

    UVector3 toleranceVector;
    toleranceVector.Set(fTolerance);

    for (int i = 0; i < numNodes; ++i)
    {
      VUSolid& solid = *solids[i];
      //UTransform3D& transform = *transforms[i];
      UTransform3D transform = transforms[i];
      UVector3 min, max;

      solid.Extent(min, max);
      if (solid.GetEntityType() == "Orb")
      {
        UOrb& orb = *(UOrb*) &solid;
        UVector3 orbToleranceVector;
        double tolerance = orb.GetRadialTolerance() / 2.0;
        orbToleranceVector.Set(tolerance);
        min -= orbToleranceVector;
        max += orbToleranceVector;
      }
      else
      {
        min -= toleranceVector;
        max += toleranceVector;
      }
      UUtils::TransformLimits(min, max, transform);
      fBoxes[i].hlen = (max - min) / 2;
      fBoxes[i].pos = transform.fTr;
    }
    fTotalCandidates = fBoxes.size();
  }

  
}
#endif // USOLIDSONLY

void UVoxelizer::BuildVoxelLimits(vector<VUFacet*>& facets)
{
  // "BuildVoxelLimits"'s aim is to store the coordinates of the origin as well as
  // the half lengths related to the bounding box of each node.
  // These quantities are stored in the array "fBoxes" (6 different values per
  // node.
  if (int numNodes = facets.size()) // Number of nodes in "multiUnion"
  {
    fBoxes.resize(numNodes); // Array which will store the half lengths
    fNPerSlice = 1 + (fBoxes.size() - 1) / (8 * sizeof(unsigned int));

    UVector3 toleranceVector;
    toleranceVector.Set(10 * fTolerance);

    for (int i = 0; i < numNodes; ++i)
    {
      VUFacet& facet = *facets[i];
      UVector3 min, max;
      UVector3 x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
      UVector3 extent;
      max.Set(facet.Extent(x), facet.Extent(y), facet.Extent(z));
      min.Set(-facet.Extent(-x), -facet.Extent(-y), -facet.Extent(-z));
      min -= toleranceVector;
      max += toleranceVector;
      UVector3 hlen = (max - min) / 2;
      fBoxes[i].hlen = hlen;
      fBoxes[i].pos = min + hlen;
    }
    fTotalCandidates = fBoxes.size();
  }
}


//______________________________________________________________________________
void UVoxelizer::DisplayVoxelLimits()
{
  // "DisplayVoxelLimits" displays the dX, dY, dZ, pX, pY and pZ for each node
  int numNodes = fBoxes.size();
  int oldprc = cout.precision();
  for (int i = 0; i < numNodes; ++i)
  {
    cout << setw(10) << setiosflags(ios::fixed) << setprecision(16) <<
         "    -> Node " << i + 1 <<  ":\n" <<
         "\t * [x,y,z] = " << fBoxes[i].hlen <<
         "\t * [x,y,z] = " << fBoxes[i].pos << "\n";
  }
  cout.precision(oldprc) ;
}

//______________________________________________________________________________
void UVoxelizer::CreateSortedBoundary(vector<double>& boundary, int axis)
{
  // "CreateBoundaries"'s aim is to determine the slices induced by the bounding fBoxes,
  // along each axis. The created boundaries are stored in the array "boundariesRaw"
  int numNodes = fBoxes.size(); // Number of nodes in structure of "UMultiUnion" type
  // Determination of the boundaries along x, y and z axis:
  for (int i = 0 ; i < numNodes; ++i)
  {
    // For each node, the boundaries are created by using the array "fBoxes"
    // built in method "BuildVoxelLimits":
    double p = fBoxes[i].pos[axis], d = fBoxes[i].hlen[axis];
    // x boundaries
    boundary[2 * i] = p - d;
    boundary[2 * i + 1] = p + d;
  }
  sort(boundary.begin(), boundary.end());
}

void UVoxelizer::BuildBoundaries()
{
  // "SortBoundaries" orders the boundaries along each axis (increasing order)
  // and also does not take into account redundant boundaries, ie if two boundaries
  // are separated by a distance strictly inferior to "tolerance".
  // The sorted boundaries are respectively stored in:
  //              * boundaries[0..2]

  // In addition, the number of elements contained in the three latter arrays are
  // precise thanks to variables: boundariesCountX, boundariesCountY and boundariesCountZ.

  if (int numNodes = fBoxes.size())
  {
    const double tolerance = fTolerance / 100.0; // Minimal distance to discriminate two boundaries.
    vector<double> sortedBoundary(2 * numNodes);

    int considered;

    for (int j = 0; j <= 2; ++j)
    {
      CreateSortedBoundary(sortedBoundary, j);
      vector<double>& boundary = fBoundaries[j];
      boundary.clear();

      considered = 0;

      for (int i = 0 ; i < 2 * numNodes; ++i)
      {
        double newBoundary = sortedBoundary[i];
        int size = boundary.size();
        if (!size || abs(boundary[size - 1] - newBoundary) > tolerance)
        {
          considered++;
          {
            boundary.push_back(newBoundary);
            continue;
          }
        }
        // If two successive boundaries are too close from each other, only the first one is considered
        {
          // cout << "Skipping boundary [" << j << "] : " << i << endl;
        }
      }

      int n = boundary.size();
      int max = 100000;
      if (n > max / 2)
      {
        int skip = n / (max / 2); // n has to be 2x bigger then 50.000. therefore only from 100.000 re
        vector<double> reduced;
        for (int i = 0; i < n; ++i)
        {
          // 50 ok for 2k, 1000, 2000
          int size = boundary.size();
          if (i % skip == 0 || i == 0 || i == size - 1) // this condition of merging boundaries was wrong, it did not count with right part, which can be completely ommited and not included in final consideration. Now should be OK
            reduced.push_back(boundary[i]);
        }
        boundary = reduced;
      }
    }
  }
}

void UVoxelizer::DisplayBoundaries()
{
  char axis[3] = {'X', 'Y', 'Z'};
  for (int i = 0; i <= 2; ++i)
  {
    cout << " * " << axis[i] << " axis:" << endl << "    | ";
    DisplayBoundaries(fBoundaries[i]);
  }
}

//______________________________________________________________________________
void UVoxelizer::DisplayBoundaries(vector<double>& boundaries)
{
  // Prints the positions of the boundaries of the slices on the three axis:

  int count = boundaries.size();
  int oldprc = cout.precision();
  for (int i = 0; i < count; ++i)
  {
    cout << setw(10) << setiosflags(ios::fixed) << setprecision(16) << boundaries[i];
    //      printf("%19.15f ", boundaries[0][i]);
    //      cout << boundaries[0][i] << " ";
    if (i != count - 1) cout << "-> ";
  }
  cout << "|" << endl << "Number of boundaries: " << count << endl;
  cout.precision(oldprc) ;
}

void UVoxelizer::BuildBitmasks(std::vector<double> boundaries[], UBits bitmasks[])
{
  // "BuildListNodes" stores in the bitmasks solids present in each slice along an axis.
  //  const double localTolerance = 0;
  int numNodes = fBoxes.size();
  int bitsPerSlice = GetBitsPerSlice();

  for (int k = 0; k < 3; k++)
  {
    int total = 0;
    vector<double>& boundary = boundaries[k];
    int voxelsCount = boundary.size() - 1;
    UBits& bitmask = bitmasks[k];

    if (bitmasks != NULL)
    {
      bitmask.Clear();
      bitmask.SetBitNumber(voxelsCount * bitsPerSlice - 1, false); // it is here so we can set the maximum number of bits. this line will rellocate the memory and set all to zero
    }
    vector<int>& candidatesCount = fCandidatesCounts[k];
    candidatesCount.resize(voxelsCount);

    for (int i = 0 ; i < voxelsCount; ++i) candidatesCount[i] = 0;

    // Loop on the nodes, number of slices per axis
    for (int j = 0 ; j < numNodes; j++)
    {
      // Determination of the minimum and maximum position along x
      // of the bounding boxe of each node:
      double p = fBoxes[j].pos[k], d = fBoxes[j].hlen[k];

      double min = p - d; // - localTolerance;
      double max = p + d; // + localTolerance;

      int i = BinarySearch(boundary, min);
      if (i < 0) i = 0;

      do
      {
        if (bitmasks != NULL)
          bitmask.SetBitNumber(i * bitsPerSlice + j);

        candidatesCount[i]++;
        total++;
        i++;
      }
      while (max > boundary[i] && i < voxelsCount);
    }

#ifdef DEBUG
    string fileName = "reduced" + UUtils::ToString(k);
    UUtils::SaveVectorToExternalFile(candidatesCount, fileName);
#endif

    /*
    if (voxelsCount < 20000)
    {
    int total2 = 0;
    for(int i = 0 ; i < voxelsCount; ++i)
    {
    candidatesCount[i] = 0;
    // Loop on the nodes, number of slices per axis
    for(int j = 0 ; j < numNodes; j++)
    {
    // Determination of the minimum and maximum position along x
    // of the bounding boxe of each node:
    double p = fBoxes[j].pos[k], d = fBoxes[j].hlen[k];
    double leftBoundary = boundary[i];
    double rightBoundary = boundary[i+1];
    double min = p - d; // - localTolerance;
    double max = p + d; // + localTolerance;
    if (min < rightBoundary && max > leftBoundary)
    {
    // Storage of the number of the node in the array:
    bitmask.SetBitNumber(i*bitsPerSlice+j);
    candidatesCount[i]++;
    total2++;
    }
    }
    }
    if (total2 != total)
    total2 = total;
    }
    */
  }
#ifdef USPECSDEBUG
  cout << "Build list nodes completed" << endl;
#endif
}

//______________________________________________________________________________
string UVoxelizer::GetCandidatesAsString(const UBits& bits)
{
  // Decodes the candidates in mask as string.
  stringstream ss;
  int numNodes = fBoxes.size();

  for (int i = 0; i < numNodes; ++i)
    if (bits.TestBitNumber(i)) ss << i + 1 << " ";
  const string& result = ss.str();
  return result;
}


//______________________________________________________________________________
void UVoxelizer::DisplayListNodes()
{
  char axis[3] = {'X', 'Y', 'Z'};
  // Prints which solids are present in the slices previously elaborated.
  int size = 8 * sizeof(int) * fNPerSlice;
  UBits bits(size);

  for (int j = 0; j <= 2; ++j)
  {
    cout << " * " << axis[j] << " axis:" << endl;
    int count = fBoundaries[j].size();
    for (int i = 0; i < count - 1; ++i)
    {
      cout << "    Slice #" << i + 1 << ": [" << fBoundaries[j][i] << " ; " << fBoundaries[j][i + 1] << "] -> ";
      bits.Set(size, (const char*)fBitmasks[j].fAllBits + i * fNPerSlice * sizeof(int));
      string result = GetCandidatesAsString(bits);
      cout << "[ " << result.c_str() << "]  " << endl;
    }
  }
}

void UVoxelizer::BuildBoundingBox()
{
  UVector3 min(fBoundaries[0].front(), fBoundaries[1].front(), fBoundaries[2].front());
  UVector3 max(fBoundaries[0].back(), fBoundaries[1].back(), fBoundaries[2].back());
  BuildBoundingBox(min, max);
}

void UVoxelizer::BuildBoundingBox(UVector3& amin, UVector3& amax, double tolerance)
{
  for (int i = 0; i <= 2; ++i)
  {
    double min = amin[i];
    double max = amax[i];
    fBoundingBoxSize[i] = (max - min) / 2 + tolerance * 0.5;
    fBoundingBoxCenter[i] = min + fBoundingBoxSize[i];
  }
  fBoundingBox = UBox("", fBoundingBoxSize.x(), fBoundingBoxSize.y(), fBoundingBoxSize.z());
}


// algorithm -

// in order to get balanced voxels, merge should always unite those regions, where the number of voxels is least
// the number

// we will keep sorted list (std::set) with all voxels. there will be comparator function between two voxels,
// which will tell if voxel is less by looking at his right neighbor.
// first, we will add all the voxels into the tree.
// we will be pick the first item in the tree, merging it, adding the right merged voxel into the a list for
// future reduction (fBitmasks will be rebuilded later, therefore they need not to be updated).
// the merged voxel need to be added to the tree again, so it's position would be updated


class UVoxelComparator
{
  public:

    vector<UVoxelInfo>& fVoxels;

    UVoxelComparator(vector<UVoxelInfo>& voxels) : fVoxels(voxels)
    {

    }

    bool operator()(int l, int r)
    {
      UVoxelInfo& lv = fVoxels[l], &rv = fVoxels[r];
      int left = lv.count +  fVoxels[lv.next].count;
      int right = rv.count + fVoxels[rv.next].count;;
      return (left == right) ? l < r : left < right;
    }
};

void UVoxelizer::SetReductionRatio(int maxVoxels, UVector3& reductionRatio)
{
  double maxTotal = (double) fCandidatesCounts[0].size() * fCandidatesCounts[1].size() * fCandidatesCounts[2].size();

  if (maxVoxels > 0 && maxVoxels < maxTotal)
  {
    double ratio = (double) maxVoxels / maxTotal;
    ratio = pow(ratio, 1. / 3.);
    if (ratio > 1) ratio = 1;
    reductionRatio.Set(ratio);
  }
}

void UVoxelizer::BuildReduceVoxels(vector<double> boundaries[], UVector3 reductionRatio)
{
  for (int k = 0; k <= 2; ++k)
  {
    vector<int>& candidatesCount = fCandidatesCounts[k];
    int max = candidatesCount.size();
    vector<UVoxelInfo> voxels(max);
    UVoxelComparator comp(voxels);
    set<int, UVoxelComparator> voxelSet(comp);
    vector<int> mergings;

    for (int j = 0; j < max; ++j)
    {
      UVoxelInfo& voxel = voxels[j];
      voxel.count = candidatesCount[j];
      voxel.previous = j - 1;
      voxel.next = j + 1;
      voxels[j] = voxel;
    }
    //    voxels[max - 1].count = 99999999999;

    for (int j = 0; j < max - 1; ++j) voxelSet.insert(j); // we go to size-1 to make sure we will not merge the last element

    double reduction = reductionRatio[k];
    if (reduction != 0)
    {
      int count = 0, currentCount;
      while ((currentCount = voxelSet.size()) > 2)
      {
        double currentRatio = 1 - (double) count / max;
        if (currentRatio <= reduction && currentCount <= 1000)
          break;
        const int pos = *voxelSet.begin();
        mergings.push_back(pos + 1);

        UVoxelInfo& voxel = voxels[pos];
        UVoxelInfo& nextVoxel = voxels[voxel.next];

        if (voxelSet.erase(pos) != 1)
        {
          ;//  k = k;
        }
        if (voxel.next != max - 1)
          if (voxelSet.erase(voxel.next) != 1)
          {
            // k = k;
          }
        if (voxel.previous != -1)
          if (voxelSet.erase(voxel.previous) != 1)
          {
            // k = k;
          }
        nextVoxel.count += voxel.count;
        voxel.count = 0;
        nextVoxel.previous = voxel.previous;

        if (voxel.next != max - 1)
          voxelSet.insert(voxel.next);

        if (voxel.previous != -1)
        {
          voxels[voxel.previous].next = voxel.next;
          voxelSet.insert(voxel.previous);
        }
        count++;
      }
    }

    // for (int i = 0; i < max; ++i) cout << voxels[i].count << ", ";

    if (mergings.size())
    {
      sort(mergings.begin(), mergings.end());

      const vector<double>& boundary = boundaries[k];
      int mergingsSize = mergings.size();
      vector<double> reducedBoundary;
      int skip = mergings[0], i = 0;
      max = boundary.size();
      for (int j = 0; j < max; ++j)
      {
        if (j != skip)
          reducedBoundary.push_back(boundary[j]);
        else if (++i < mergingsSize)
        {
          skip = mergings[i];
        }
      }
      boundaries[k] = reducedBoundary;
    }
  }
}


void UVoxelizer::BuildReduceVoxels2(vector<double> boundaries[], UVector3 reductionRatio)
{
  for (int k = 0; k <= 2; ++k)
  {
    vector<int>& candidatesCount = fCandidatesCounts[k];
    int max = candidatesCount.size();
    int total = 0;
    for (int i = 0; i < max; ++i) total += candidatesCount[i];

    double reduction = reductionRatio[k];
    if (reduction == 0)
      break;

    int destination = (int)(reduction * max) + 1;
    if (destination > 1000) destination = 1000;
    if (destination < 2) destination = 2;
    double average = ((double)total / max) / reduction;

    vector<int> mergings;

    vector<double>& boundary = boundaries[k];
    vector<double> reducedBoundary(destination);

    int sum = 0, cur = 0;
    for (int i = 0; i < max; ++i)

    {
      sum += candidatesCount[i];
      if (sum > average * (cur + 1) || i == 0)
      {
        double val = boundary[i];
        reducedBoundary[cur] = val;
        cur++;
        if (cur == destination)
          break;
      }
    }
    reducedBoundary[destination - 1] = boundary[max];
    boundaries[k] = reducedBoundary;
  }
}

#ifdef USOLIDSONLY
//______________________________________________________________________________
void UVoxelizer::Voxelize(vector<VUSolid*>& solids, vector<UTransform3D>& transforms)
{
  
  BuildVoxelLimits(solids, transforms);

  BuildBoundaries();

  BuildBitmasks(fBoundaries, fBitmasks);

  BuildBoundingBox();

  BuildEmpty(); // these does not work well for multi-union, actually only makes performance slower, these are only pre-calculated but not used by multi-union

//  fBoxes.resize(0); // for multiunion, fBoxes are actually needed

  for (int i = 0; i < 3; ++i) fCandidatesCounts[i].resize(0);
}
#endif // USOLIDSONLY

void UVoxelizer::CreateMiniVoxels(std::vector<double> boundaries[], UBits bitmasks[])
{
  vector<int> voxel(3), maxVoxels(3);
  for (int i = 0; i <= 2; ++i) maxVoxels[i] = boundaries[i].size();

  UVector3 point;
  for (voxel[2] = 0; voxel[2] < maxVoxels[2] - 1; ++voxel[2])
  {
    for (voxel[1] = 0; voxel[1] < maxVoxels[1] - 1; ++voxel[1])
    {
      for (voxel[0] = 0; voxel[0] < maxVoxels[0] - 1; ++voxel[0])
      {
        vector<int> candidates;
        if (GetCandidatesVoxelArray(voxel, bitmasks, candidates, NULL))
        {
          // find a box for corresponding non-empty voxel
          UVoxelBox box;
          for (int i = 0; i <= 2; ++i)
          {
            int index = voxel[i];
            const vector<double>& boundary = boundaries[i];
            double hlen = 0.5 * (boundary[index + 1] - boundary[index]);
            box.hlen[i] = hlen;
            box.pos[i] = boundary[index] + hlen;
          }
          fVoxelBoxes.push_back(box);
          vector<int>(candidates).swap(candidates);
          fVoxelBoxesCandidates.push_back(candidates);
        }
      }
    }
  }
}

void UVoxelizer::Voxelize(vector<VUFacet*>& facets)
{
  int maxVoxels = fMaxVoxels;
  UVector3 reductionRatio = fReductionRatio;

  int size = facets.size();
  if (size < 10)
    for (int i = 0; i < (int) facets.size(); ++i)
      if (facets[i]->GetNumberOfVertices() > 3) size++;

  if ((size >= 10 || maxVoxels > 0) && maxVoxels != 0 && maxVoxels != 1)
  {
    BuildVoxelLimits(facets);

    BuildBoundaries();

    BuildBitmasks(fBoundaries, NULL);

    if (maxVoxels < 0 && reductionRatio == UVector3())
    {
      /*
      if (fTotalCandidates < 1000)
        maxVoxels = 10 * fTotalCandidates;
      else
      {
        if (fTotalCandidates < 10000) maxVoxels = 10000;
        else maxVoxels = fTotalCandidates;
      }
      if (fTotalCandidates > 1000000) maxVoxels = 1000000;
      */

      maxVoxels = fTotalCandidates;
      if (fTotalCandidates > 1000000) maxVoxels = 1000000;
    }

    SetReductionRatio(maxVoxels, reductionRatio);

    fCountOfVoxels = CountVoxels(fBoundaries);

#ifdef USPECSDEBUG
    cout << "Total number of voxels: " << fCountOfVoxels << endl;
#endif

    BuildReduceVoxels2(fBoundaries, reductionRatio);

    fCountOfVoxels = CountVoxels(fBoundaries);

#ifdef USPECSDEBUG
    cout << "Total number of voxels after reduction: " << fCountOfVoxels << endl;
#endif

    BuildBitmasks(fBoundaries, fBitmasks);

    UVector3 reductionRatioMini;

    UBits bitmasksMini[3];

    // section for building mini voxels
    std::vector<double> miniBoundaries[3];

    for (int i = 0; i <= 2; ++i) miniBoundaries[i] = fBoundaries[i];

    int voxelsCountMini = (fCountOfVoxels >= 1000) ? 100 : fCountOfVoxels / 10;
    //  if (voxelCountMini < 8) voxelCountMini = 8;

//    voxelsCountMini = 1;

    SetReductionRatio(voxelsCountMini, reductionRatioMini);

    BuildReduceVoxels(miniBoundaries, reductionRatioMini);

#ifdef USPECSDEBUG
    int total = CountVoxels(miniBoundaries);
    cout << "Total number of mini voxels: " << total << endl;
#endif

    BuildBitmasks(miniBoundaries, bitmasksMini);

    CreateMiniVoxels(miniBoundaries, bitmasksMini);

    BuildEmpty();

    // deallocate fields unnecessary during runtime
    fBoxes.resize(0);

    for (int i = 0; i < 3; ++i)
    {
      fCandidatesCounts[i].resize(0);
      fBitmasks[i].Clear();
    }
  }
}


//______________________________________________________________________________
// "GetCandidates" should compute which solids are possibly contained in
// the voxel defined by the three slices characterized by the passed indexes.
void UVoxelizer::GetCandidatesVoxel(vector<int>& voxels)
{
  cout << "   Candidates in voxel [" << voxels[0] << " ; " << voxels[1] << " ; " << voxels[2] << "]: ";
  vector<int> candidates;
  int count = GetCandidatesVoxelArray(voxels, candidates);
  cout << "[ ";
  for (int i = 0; i < count; ++i) cout << candidates[i];
  cout << "]  " << endl;
}

#ifdef USOLIDSONLY

inline void findComponents2(unsigned int mask, vector<int>& list, int i)
{
  for (int bit = 0; bit < (int)(8 * sizeof(unsigned int)); bit++)
  {
    if (mask & 1)
      list.push_back(8 * sizeof(unsigned int)*i + bit);
    if (!(mask >>= 1)) break; // new
  }
}

inline void findComponents3(unsigned int mask, vector<int>& list, int i)
{
  // nejrychlejsi asi bude: ve for cyclu traversovat pres vsechny byty:
  // 1. voendovanou hodnotu pres 0xFF, aby zustal jen posledni byte
  // 2. pokud 0, continue
  // 3. najit prvni nastaveny bit pres lookup table
  // 4. pricist 8*j k nalezenemu bitu, ulozit do seznamu
  // 5. odecist 1 << bit, dokud neni *-*nula pokracovat na 3
  // 6. posunout hodnotu do prava >> 8, pokracovat na 1

  for (int byte = 0; byte < (int)(sizeof(unsigned int)); byte++)
  {
    if (int maskByte = mask & 0xFF)
    {
      do
      {
        static const int firstBits[256] =
        {
          8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
          4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
        };

        int bit = firstBits[maskByte];
        list.push_back(8 * (sizeof(unsigned int)*i + byte) + bit);
        maskByte -= 1 << bit;
      }
      while (maskByte);
    }
    mask >>= 8;
  }
}

#endif // USOLIDSONLY

void UVoxelizer::FindComponentsFastest(unsigned int mask, vector<int>& list, int i)
{
  for (int byte = 0; byte < (int)(sizeof(unsigned int)); byte++)
  {
    if (int maskByte = mask & 0xFF)
    {
      for (int bit = 0; bit < 8; bit++)
      {
        if (maskByte & 1)
          list.push_back(8 * (sizeof(unsigned int)*i + byte) + bit);
        if (!(maskByte >>= 1)) break;
      }
    }
    mask >>= 8;
  }
}


//______________________________________________________________________________
// Method returning the candidates corresponding to the passed point
int UVoxelizer::GetCandidatesVoxelArray(const UVector3& point, vector<int>& list, UBits* crossed) const
{
  list.clear();

  for (int i = 0; i <= 2; ++i)
    if (point[i] < fBoundaries[i].front() || point[i] >= fBoundaries[i].back())
      return 0;

  if (fTotalCandidates == 1)
  {
    list.push_back(0);
    return 1;
  }
  else
  {
    if (fNPerSlice == 1)
    {
      unsigned int mask = 0xFFffFFff;
      int slice;
      if (fBoundaries[0].size() > 2)
      {
        slice = BinarySearch(fBoundaries[0], point.x());
        if (!(mask = ((unsigned int*) fBitmasks[0].fAllBits)[slice])) return 0;
      }
      if (fBoundaries[1].size() > 2)
      {
        slice = BinarySearch(fBoundaries[1], point.y());
        if (!(mask &= ((unsigned int*) fBitmasks[1].fAllBits)[slice]
             )) return 0;
      }
      if (fBoundaries[2].size() > 2)
      {
        slice = BinarySearch(fBoundaries[2], point.z());
        if (!(mask &= ((unsigned int*) fBitmasks[2].fAllBits)[slice]
             )) return 0;
      }
      if (crossed && (!(mask &= ~((unsigned int*)crossed->fAllBits)[0]))) return 0;

      FindComponentsFastest(mask, list, 0);
    }
    else
    {
      unsigned int* masks[3], mask; // masks for X,Y,Z axis
      for (int i = 0; i <= 2; ++i)
      {
        int slice = BinarySearch(fBoundaries[i], point[i]);
        //      if (slice < 0 || slice == fBoundaries[i].size()-1) return 0; // not neccesary anymore
        masks[i] = ((unsigned int*) fBitmasks[i].fAllBits) + slice * fNPerSlice;
      }
      unsigned int* maskCrossed = crossed ? (unsigned int*)crossed->fAllBits : NULL;

      for (int i = 0 ; i < fNPerSlice; ++i)
      {
        // Logic "and" of the masks along the 3 axes x, y, z:
        // removing "if (!" and ") continue" => slightly slower
        if (!(mask = masks[0][i])) continue;
        if (!(mask &= masks[1][i])) continue;
        if (!(mask &= masks[2][i])) continue;
        if (maskCrossed && !(mask &= ~maskCrossed[i])) continue;

        FindComponentsFastest(mask, list, i);
      }
    }
  }
  return list.size();
}





int UVoxelizer::GetCandidatesVoxelArray(const vector<int>& voxels, const UBits bitmasks[], vector<int>& list, UBits* crossed) const
{
  list.clear();

  if (fTotalCandidates == 1)
  {
    list.push_back(0);
    return 1;
  }
  else
  {
    if (fNPerSlice == 1)
    {
      unsigned int mask;
      if (!(mask = ((unsigned int*) bitmasks[0].fAllBits)[voxels[0]]
           )) return 0;
      if (!(mask &= ((unsigned int*) bitmasks[1].fAllBits)[voxels[1]]
           )) return 0;
      if (!(mask &= ((unsigned int*) bitmasks[2].fAllBits)[voxels[2]]
           )) return 0;
      if (crossed && (!(mask &= ~((unsigned int*)crossed->fAllBits)[0]))) return 0;

      FindComponentsFastest(mask, list, 0);
    }
    else
    {
      unsigned int* masks[3], mask; // masks for X,Y,Z axis
      for (int i = 0; i <= 2; ++i)
        masks[i] = ((unsigned int*) bitmasks[i].fAllBits) + voxels[i] * fNPerSlice;

      unsigned int* maskCrossed = crossed ? (unsigned int*)crossed->fAllBits : NULL;

      for (int i = 0 ; i < fNPerSlice; ++i)
      {
        // Logic "and" of the masks along the 3 axes x, y, z:
        // removing "if (!" and ") continue" => slightly slower
        if (!(mask = masks[0][i])) continue;
        if (!(mask &= masks[1][i])) continue;
        if (!(mask &= masks[2][i])) continue;
        if (maskCrossed && !(mask &= ~maskCrossed[i])) continue;

        FindComponentsFastest(mask, list, i);
      }
    }
  }
  return list.size();
}



//______________________________________________________________________________
// Method returning the candidates corresponding to the passed point
int UVoxelizer::GetCandidatesVoxelArray(const vector<int>& voxels, vector<int>& list, UBits* crossed) const
{
  return GetCandidatesVoxelArray(voxels, fBitmasks, list, crossed);
}



bool UVoxelizer::Contains(const UVector3& point) const
{
  for (int i = 0; i < 3; ++i)
    if (point[i] < fBoundaries[i].front() || point[i] > fBoundaries[i].back())
      return false;
  return true;
}



double UVoxelizer::DistanceToFirst(const UVector3& point, const UVector3& direction) const
{
  UVector3 pointShifted = point - fBoundingBoxCenter;
  double shift = fBoundingBox.DistanceToIn(pointShifted, direction);
  return shift;
}

double UVoxelizer::SafetyToBoundingBox(const UVector3& point) const
{
  UVector3 pointShifted = point - fBoundingBoxCenter;
  double shift = MinDistanceToBox(pointShifted, fBoundingBoxSize);
  return shift;
}

double UVoxelizer::MinDistanceToBox(const UVector3& aPoint, const UVector3& f)
{
  // Estimates the isotropic safety from a point outside the current solid to any
  // of its surfaces. The algorithm may be accurate or should provide a fast
  // underestimate.
  double safe, safx, safy, safz;
  safe = safx = -f.x() + std::abs(aPoint.x());
  safy = -f.y() + std::abs(aPoint.y());
  if (safy > safe) safe = safy;
  safz = -f.z() + std::abs(aPoint.z());
  if (safz > safe) safe = safz;
  if (safe < 0.0) return 0.0; // point is inside
//  if (!aAccurate) return safe;
  double safsq = 0.0;
  int count = 0;
  if (safx > 0)
  {
    safsq += safx * safx;
    count++;
  }
  if (safy > 0)
  {
    safsq += safy * safy;
    count++;
  }
  if (safz > 0)
  {
    safsq += safz * safz;
    count++;
  }
  if (count == 1) return safe;
  return std::sqrt(safsq);
}


double UVoxelizer::DistanceToNext(const UVector3& point, const UVector3& direction, vector<int>& curVoxel) const
{
  double shift = UUtils::kInfinity;
  int cur = 0; // the smallest index, which would be than increased

  for (int i = 0; i <= 2; ++i)
  {
    // Looking for the next voxels on the considered direction X,Y,Z axis
    const vector<double>& boundary = fBoundaries[i];
    int index = curVoxel[i];
    if (direction[i] >= 1e-10)
    {
      ++index;
//      if (boundary[++index] - point[i] < fTolerance)

//        if (++index >= (int) boundary.size() - 1)
//          exit = true;

// make sure shift would be non-zero => not needed anymore

//        if (++index >= (int) boundary.size()) // the shift is not important, only is important to increase the boundary index, even if it would be zero, than we would in next steps always increase another dimension, so it would continue alright
//          continue;
    }
    else
    {
      if (direction[i] <= -1e-10)
      {
//        if (point[i] - boundary[index] < fTolerance) // make sure shift would be non-zero
//          if (--index < 0)
//            continue;

//          if (--index < 0)
//            exit = true;
      }
      else
        continue;
    }
    double dif = boundary[index] - point[i];
    double distance = dif / direction[i];

    if (shift > distance)
    {
      shift = distance;
      cur = i;
    }
  }

  if (shift != UUtils::kInfinity)
  {
    // updating current voxel using the index corresponding to the closest voxel boundary on the ray
    if (direction[cur] > 0)
    {
      if (++curVoxel[cur] >= (int) fBoundaries[cur].size() - 1)
        shift = UUtils::kInfinity;
    }
    else
    {
      if (--curVoxel[cur] < 0)
        shift = UUtils::kInfinity;
    }
  }

  return shift;
}




void UVoxelizer::SetMaxVoxels(int max)
{
  fMaxVoxels = max;
  fReductionRatio.Set(0);
}

void UVoxelizer::SetMaxVoxels(const UVector3& ratioOfReduction)
{
  fMaxVoxels = -1;
  fReductionRatio = ratioOfReduction;
}

int UVoxelizer::AllocatedMemory()
{
  int size = fEmpty.GetNbytes();
  size += fBoxes.capacity() * sizeof(UVoxelBox);
  size += sizeof(double) * (fBoundaries[0].capacity() + fBoundaries[1].capacity() + fBoundaries[2].capacity());
  size += sizeof(int) * (fCandidatesCounts[0].capacity() + fCandidatesCounts[1].capacity() + fCandidatesCounts[2].capacity());
  size += fBitmasks[0].GetNbytes() + fBitmasks[1].GetNbytes() + fBitmasks[2].GetNbytes();

  int csize = fCandidates.size();
  for (int i = 0; i < csize; ++i)
    size += sizeof(vector<int>) + fCandidates[i].capacity() * sizeof(int);

  return size;
}

void UVoxelizer::SetDefaultVoxelsCount(int count)
{
  fDefaultVoxelsCount = count;
}

int UVoxelizer::GetDefaultVoxelsCount()
{
  return fDefaultVoxelsCount;
}

int UVoxelizer::fDefaultVoxelsCount = -1;

#ifdef USOLIDSONLY

/*
// THIS METHOD IS NOT NECCESSARY. destructors are called automatically. it is enough to call destructor of
// object which contains the voxel finder
//
inline void DeleteObjects()
{
for (int k = 0; k < 3; k++)
{
std::vector<double> &boundary = fBoundaries[k];
boundary.resize(0);
UBits &bitmask = fBitmasks[k];
bitmask.Clear();
bitmask.SetBitNumber(1, false); // it is here so we can set the maximum
}
}
*/

/*
UVoxelCandidatesIterator::UVoxelCandidatesIterator(const UVoxelizer &f, const UVector3 &point) : curInt(-1), curBit((int) (8*sizeof(unsigned int))), nextAvailable(true)
{
  carNodes = f.fBoxes.size();
  if (carNodes > 1)
  {
    n = 1+(carNodes-1)/(8*sizeof(unsigned int));
    int sliceX = UVoxelizer::BinarySearch(f.fBoundaries[0], point.x());
    int sliceY = UVoxelizer::BinarySearch(f.fBoundaries[1], point.y());
    int sliceZ = UVoxelizer::BinarySearch(f.fBoundaries[2], point.z());

    unsigned int *maskX = ((unsigned int *) f.fBitmasks[0].fAllBits) + sliceX*n;
    maskXLeft = (sliceX && point.x() == f.fBoundaries[0][sliceX]) ? maskX - n : NULL;

    unsigned int *maskY = ((unsigned int *) f.fBitmasks[1].fAllBits) + sliceY*n;
    maskYLeft = (sliceY && point.y() == f.fBoundaries[1][sliceY]) ? maskY - n : NULL;
    unsigned int *maskZ = ((unsigned int *) f.fBitmasks[2].fAllBits) + sliceZ*n;
    maskZLeft = (sliceZ && point.z() == f.fBoundaries[2][sliceZ]) ? maskZ - n  : NULL;

    if (sliceX == -1 || sliceX == f.fBoundaries[0].back()) nextAvailable = false;
    if (sliceY == -1 || sliceY == f.fBoundaries[1].back()) nextAvailable = false;
    if (sliceZ == -1 || sliceZ == f.fBoundaries[2].back()) nextAvailable = false;
  }
  else
  {
    if (point.x() < f.fBoundaries[0].front() || point.x() > f.fBoundaries[0].back()) nextAvailable = false;
    if (point.y() < f.fBoundaries[1].front() || point.y() > f.fBoundaries[1].back()) nextAvailable = false;
    if (point.z() < f.fBoundaries[2].front() || point.z() > f.fBoundaries[2].back()) nextAvailable = false;
  }
}

int UVoxelCandidatesIterator::Next()
{
  if (!nextAvailable) return -1;

  if (carNodes == 1)
  {
    nextAvailable = false;
    return 0;
  }
  else
  {
    do
    {
      while (curBit == (int) (8*sizeof(unsigned int)))
      {
        if (curInt++ >= n)
        {
          nextAvailable = false;
          return -1;
        }
        // Logic "and" of the masks along the 3 axes x, y, z:
        // removing "if (!" and ") continue" => slightly slower
        if (!(mask = maskXLeft ? maskX[curInt] | maskXLeft[curInt] : maskX[curInt]) ||
          !(mask &= maskYLeft ? maskY[curInt] | maskYLeft[curInt] : maskY[curInt]) ||
          !(mask &= maskZLeft ? maskZ[curInt] | maskZLeft[curInt] : maskZ[curInt]))
        {
          continue;
        }
        curBit = 0;
      }
      int shifted = 1 << curBit; // new
      if (mask & shifted)
      {
        int res = 8*sizeof(unsigned int)*curInt+curBit;
        if (!(mask -= shifted))
          curBit = (int) (8*sizeof(unsigned int));

        return res;
      }
      curBit++;
    }
    while (nextAvailable);
    return -1;
  }
}
*/

#endif // USOLIDSONLY
