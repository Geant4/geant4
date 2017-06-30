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
// G4Voxelizer implementation
//
// History:
// 19.10.12 Marek Gayer, created
// --------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <set>

#include "G4VSolid.hh" 

#include "G4Orb.hh"
#include "G4Voxelizer.hh"
#include "G4SolidStore.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeometryTolerance.hh"
#include "G4CSGSolid.hh"
#include "G4Orb.hh"
#include "G4Types.hh"
#include "geomdefs.hh"

using namespace std;

G4ThreadLocal G4int G4Voxelizer::fDefaultVoxelsCount = -1;

//______________________________________________________________________________
G4Voxelizer::G4Voxelizer()
  : fBoundingBox("VoxBBox", 1, 1, 1)
{
  fCountOfVoxels = fNPerSlice = fTotalCandidates = 0;
  
  fTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  SetMaxVoxels(fDefaultVoxelsCount); 

  G4SolidStore::GetInstance()->DeRegister(&fBoundingBox);
}

//______________________________________________________________________________
G4Voxelizer::~G4Voxelizer()
{
}

//______________________________________________________________________________
void G4Voxelizer::BuildEmpty()
{
  // by reserving the size of candidates, we would avoid reallocation of 
  // the vector which could cause fragmentation
  //
  std::vector<G4int> xyz(3), max(3), candidates(fTotalCandidates);
  const std::vector<G4int> empty(0);

  for (G4int i = 0; i <= 2; i++) max[i] = fBoundaries[i].size();
  unsigned int size = max[0] * max[1] * max[2];

  fEmpty.Clear();
  fEmpty.ResetBitNumber(size-1);
  fEmpty.ResetAllBits(true);

  for (xyz[2] = 0; xyz[2] < max[2]; ++xyz[2])
  {
    for (xyz[1] = 0; xyz[1] < max[1]; ++xyz[1])
    {
      for (xyz[0] = 0; xyz[0] < max[0]; ++xyz[0])
      {
        if (GetCandidatesVoxelArray(xyz, candidates)) 
        {
          G4int index = GetVoxelsIndex(xyz);
          fEmpty.SetBitNumber(index, false);

          // rather than assigning directly with:
          // "fCandidates[index] = candidates;", in an effort to ensure that 
          // capacity would be just exact, we rather use following 3 lines
          //
          std::vector<G4int> &c = (fCandidates[index] = empty);
          c.reserve(candidates.size());
          c.assign(candidates.begin(), candidates.end());
        }
      }
    }
  }
#ifdef G4SPECSDEBUG
  G4cout << "Non-empty voxels count: " << fCandidates.size() << G4endl;
#endif
}

//______________________________________________________________________________
void G4Voxelizer::BuildVoxelLimits(std::vector<G4VSolid*>& solids,
                                   std::vector<G4Transform3D>& transforms)
{
  G4Rotate3D rot;
  G4Translate3D transl ;
  G4Scale3D scale;

  // "BuildVoxelLimits"'s aim is to store the coordinates of the origin as
  // well as the half lengths related to the bounding box of each node.
  // These quantities are stored in the array "fBoxes" (6 different values per
  // node
  //
  if (G4int numNodes = solids.size()) // Number of nodes in "multiUnion"
  {
    fBoxes.resize(numNodes); // Array which will store the half lengths
    fNPerSlice = 1 + (fBoxes.size() - 1) / (8 * sizeof(unsigned int));

    // related to a particular node, but also
    // the coordinates of its origin

    G4ThreeVector toleranceVector(fTolerance,fTolerance,fTolerance);

    for (G4int i = 0; i < numNodes; ++i)
    {
      G4VSolid& solid = *solids[i];
      G4Transform3D transform = transforms[i];
      G4ThreeVector min, max;

      solid.BoundingLimits(min, max);
      if (solid.GetEntityType() == "G4Orb")
      {
        G4Orb& orb = *(G4Orb*) &solid;
        G4ThreeVector orbToleranceVector;
        G4double tolerance = orb.GetRadialTolerance() / 2.0;
        orbToleranceVector.set(tolerance,tolerance,tolerance);
        min -= orbToleranceVector;
        max += orbToleranceVector;
      }
      else
      {
        min -= toleranceVector;
        max += toleranceVector;
      }
      TransformLimits(min, max, transform);
      fBoxes[i].hlen = (max - min) / 2;
      transform.getDecomposition(scale,rot,transl); 
      fBoxes[i].pos = transl.getTranslation();
    }
    fTotalCandidates = fBoxes.size();
  }
}

//______________________________________________________________________________
void G4Voxelizer::BuildVoxelLimits(std::vector<G4VFacet *> &facets)
{
  // "BuildVoxelLimits"'s aim is to store the coordinates of the origin as well
  // as the half lengths related to the bounding box of each node.
  // These quantities are stored in the array "fBoxes" (6 different values per
  // node.

  if (G4int numNodes = facets.size()) // Number of nodes
  {
    fBoxes.resize(numNodes); // Array which will store the half lengths
    fNPerSlice = 1+(fBoxes.size()-1)/(8*sizeof(unsigned int));

    G4ThreeVector toleranceVector(10*fTolerance, 10*fTolerance, 10*fTolerance);

    for (G4int i = 0; i < numNodes; ++i)
    {
      G4VFacet &facet = *facets[i];
      G4ThreeVector min, max;
      G4ThreeVector x(1,0,0), y(0,1,0), z(0,0,1);
      G4ThreeVector extent;
      max.set (facet.Extent(x), facet.Extent(y), facet.Extent(z));
      min.set (-facet.Extent(-x), -facet.Extent(-y), -facet.Extent(-z));
      min -= toleranceVector;
      max += toleranceVector;
      G4ThreeVector hlen = (max - min) / 2;
      fBoxes[i].hlen = hlen;
      fBoxes[i].pos = min + hlen;
    }
    fTotalCandidates = fBoxes.size();
  }
}

//______________________________________________________________________________
void G4Voxelizer::DisplayVoxelLimits() const
{
  // "DisplayVoxelLimits" displays the dX, dY, dZ, pX, pY and pZ for each node

  G4int numNodes = fBoxes.size();
  G4int oldprec = G4cout.precision(16);
  for(G4int i = 0; i < numNodes; ++i)
  {
    G4cout << setw(10) << setiosflags(ios::fixed) <<
      "    -> Node " << i+1 <<  ":\n" << 
      "\t * [x,y,z] = " << fBoxes[i].hlen <<
      "\t * [x,y,z] = " << fBoxes[i].pos << "\n";
  }
  G4cout.precision(oldprec);
}

//______________________________________________________________________________
void G4Voxelizer::CreateSortedBoundary(std::vector<G4double> &boundary,
                                       G4int axis)
{
  // "CreateBoundaries"'s aim is to determine the slices induced by the
  // bounding fBoxes, along each axis. The created boundaries are stored
  // in the array "boundariesRaw"

  G4int numNodes = fBoxes.size(); // Number of nodes in structure

  // Determination of the boundaries along x, y and z axis
  //
  for(G4int i = 0 ; i < numNodes; ++i)   
  {
    // For each node, the boundaries are created by using the array "fBoxes"
    // built in method "BuildVoxelLimits"
    //
    G4double p = fBoxes[i].pos[axis], d = fBoxes[i].hlen[axis];

    // x boundaries
    //
#ifdef G4SPECSDEBUG
    G4cout << "Boundary " << p - d << " - " << p + d << G4endl;
#endif
    boundary[2*i] = p - d;
    boundary[2*i+1] = p + d;
  }
  std::sort(boundary.begin(), boundary.end());
}

//______________________________________________________________________________
void G4Voxelizer::BuildBoundaries()
{
  // "SortBoundaries" orders the boundaries along each axis (increasing order)
  // and also does not take into account redundant boundaries, i.e. if two
  // boundaries are separated by a distance strictly inferior to "tolerance".
  // The sorted boundaries are respectively stored in:
  //              * boundaries[0..2]

  // In addition, the number of elements contained in the three latter arrays
  // are precise thanks to variables: boundariesCountX, boundariesCountY and
  // boundariesCountZ.

  if (G4int numNodes = fBoxes.size())
  {
    const G4double tolerance = fTolerance / 100.0;
      // Minimal distance to discriminate two boundaries.

    std::vector<G4double> sortedBoundary(2*numNodes);

    G4int considered;

    for (G4int j = 0; j <= 2; ++j)
    {
      CreateSortedBoundary(sortedBoundary, j);
      std::vector<G4double> &boundary = fBoundaries[j];
      boundary.clear();

      considered = 0;

      for(G4int i = 0 ; i < 2*numNodes; ++i)
      {
        G4double newBoundary = sortedBoundary[i];
#ifdef G4SPECSDEBUG	
        if (j == 0) G4cout << "Examining " << newBoundary << "..." << G4endl;
#endif
        G4int size = boundary.size();
        if(!size || std::abs(boundary[size-1] - newBoundary) > tolerance)
        {
          considered++;
          {
#ifdef G4SPECSDEBUG	    
            if (j == 0) G4cout << "Adding boundary " << newBoundary << "..."
                               << G4endl;
#endif	  
            boundary.push_back(newBoundary);
            continue;
          }
        }
        // If two successive boundaries are too close from each other,
        // only the first one is considered   
      }

      G4int n = boundary.size();
      G4int max = 100000;
      if (n > max/2)
      {
        G4int skip = n / (max /2); // n has to be 2x bigger then 50.000.
                                   // therefore only from 100.000 reduced
        std::vector<G4double> reduced;
        for (G4int i = 0; i < n; ++i)
        {
          // 50 ok for 2k, 1000, 2000
          G4int size = boundary.size();
          if (i % skip == 0 || i == 0 || i == size - 1)
          {
            // this condition of merging boundaries was wrong,
            // it did not count with right part, which can be
            // completely ommited and not included in final consideration.
            // Now should be OK
            //
            reduced.push_back(boundary[i]);
          }
        }
        boundary = reduced;
      }
    }
  }
}

//______________________________________________________________________________
void G4Voxelizer::DisplayBoundaries()
{
  char axis[3] = {'X', 'Y', 'Z'};
  for (G4int i = 0; i <= 2; ++i)
  {
    G4cout << " * " << axis[i] << " axis:" << G4endl << "    | ";
    DisplayBoundaries(fBoundaries[i]);
  }
}

//______________________________________________________________________________
void G4Voxelizer::DisplayBoundaries(std::vector<G4double> &boundaries)
{
  // Prints the positions of the boundaries of the slices on the three axes

  G4int count = boundaries.size();
  G4int oldprec = G4cout.precision(16);
  for(G4int i = 0; i < count; ++i)
  {
    G4cout << setw(10) << setiosflags(ios::fixed) << boundaries[i];
    if(i != count-1) G4cout << "-> ";
  }
  G4cout << "|" << G4endl << "Number of boundaries: " << count << G4endl;
  G4cout.precision(oldprec);
}

//______________________________________________________________________________
void G4Voxelizer::BuildBitmasks(std::vector<G4double> boundaries[],
                                G4SurfBits bitmasks[], G4bool countsOnly)
{
  // "BuildListNodes" stores in the bitmasks solids present in each slice
  // along an axis.

  G4int numNodes = fBoxes.size();
  G4int bitsPerSlice = GetBitsPerSlice();

  for (G4int k = 0; k < 3; ++k)
  {
    G4int total = 0;
    std::vector<G4double> &boundary = boundaries[k];
    G4int voxelsCount = boundary.size() - 1;
    G4SurfBits &bitmask = bitmasks[k];

    if (!countsOnly)
    {
      bitmask.Clear();
#ifdef G4SPECSDEBUG
      G4cout << "Allocating bitmask..." << G4endl;
#endif
      bitmask.SetBitNumber(voxelsCount*bitsPerSlice-1, false);
        // it is here so we can set the maximum number of bits. this line
        // will rellocate the memory and set all to zero
    }
    std::vector<G4int> &candidatesCount = fCandidatesCounts[k];
    candidatesCount.resize(voxelsCount);

    for(G4int i = 0 ; i < voxelsCount; ++i) { candidatesCount[i] = 0; }

    // Loop on the nodes, number of slices per axis
    //
    for(G4int j = 0 ; j < numNodes; ++j)
    {
      // Determination of the minimum and maximum position along x
      // of the bounding boxe of each node
      //
      G4double p = fBoxes[j].pos[k], d = fBoxes[j].hlen[k];

      G4double min = p - d; // - localTolerance;
      G4double max = p + d; // + localTolerance;

      G4int i = BinarySearch(boundary, min);
      if (i < 0)  { i = 0; }

      do    // Loop checking, 13.08.2015, G.Cosmo
      {
        if (!countsOnly)
        {
          bitmask.SetBitNumber(i*bitsPerSlice+j);
        }
        candidatesCount[i]++;
        total++;
        i++;
      }
      while (max > boundary[i] && i < voxelsCount);
    }
  }
#ifdef G4SPECSDEBUG
  G4cout << "Build list nodes completed." << G4endl;
#endif
}

//______________________________________________________________________________
G4String G4Voxelizer::GetCandidatesAsString(const G4SurfBits &bits) const
{
  // Decodes the candidates in mask as G4String.

  stringstream ss;
  G4int numNodes = fBoxes.size();

  for(G4int i=0; i<numNodes; ++i)
  {
    if (bits.TestBitNumber(i))  { ss << i+1 << " "; }
  }
  return ss.str();
}

//______________________________________________________________________________
void G4Voxelizer::DisplayListNodes() const
{
  // Prints which solids are present in the slices previously elaborated.

  char axis[3] = {'X', 'Y', 'Z'};
  G4int size=8*sizeof(G4int)*fNPerSlice;
  G4SurfBits bits(size);

  for (G4int j = 0; j <= 2; ++j)
  {
    G4cout << " * " << axis[j] << " axis:" << G4endl;
    G4int count = fBoundaries[j].size();
    for(G4int i=0; i < count-1; ++i)
    {
      G4cout << "    Slice #" << i+1 << ": [" << fBoundaries[j][i]
             << " ; " << fBoundaries[j][i+1] << "] -> ";
      bits.set(size,(const char *)fBitmasks[j].fAllBits+i
                                 *fNPerSlice*sizeof(G4int));
      G4String result = GetCandidatesAsString(bits);
      G4cout << "[ " << result.c_str() << "]  " << G4endl;
    }
  }
}

//______________________________________________________________________________
void G4Voxelizer::BuildBoundingBox()
{
  G4ThreeVector min(fBoundaries[0].front(),
                    fBoundaries[1].front(),
                    fBoundaries[2].front());
  G4ThreeVector max(fBoundaries[0].back(),
                    fBoundaries[1].back(),
                    fBoundaries[2].back());
  BuildBoundingBox(min, max);
}

//______________________________________________________________________________
void G4Voxelizer::BuildBoundingBox(G4ThreeVector& amin,
                                   G4ThreeVector& amax,
                                   G4double tolerance)
{
  for (G4int i = 0; i <= 2; ++i)
  {
    G4double min = amin[i];
    G4double max = amax[i];
    fBoundingBoxSize[i] = (max - min) / 2 + tolerance * 0.5;
    fBoundingBoxCenter[i] = min + fBoundingBoxSize[i];
  }
  fBoundingBox = G4Box("VoxBBox", fBoundingBoxSize.x(),
                       fBoundingBoxSize.y(), fBoundingBoxSize.z());
}

// algorithm - 

// in order to get balanced voxels, merge should always unite those regions,
// where the number of voxels is least the number.
// We will keep sorted list (std::set) with all voxels. there will be
// comparator function between two voxels, which will tell if voxel is less
// by looking at his right neighbor.
// First, we will add all the voxels into the tree.
// We will be pick the first item in the tree, merging it, adding the right
// merged voxel into the a list for future reduction (fBitmasks will be
// rebuilded later, therefore they need not to be updated).
// The merged voxel need to be added to the tree again, so it's position
// would be updated.

//______________________________________________________________________________
void G4Voxelizer::SetReductionRatio(G4int maxVoxels,
                                    G4ThreeVector &reductionRatio)
{
  G4double maxTotal = (G4double) fCandidatesCounts[0].size()
                    * fCandidatesCounts[1].size() * fCandidatesCounts[2].size();

  if (maxVoxels > 0 && maxVoxels < maxTotal)
  {
    G4double ratio = (G4double) maxVoxels / maxTotal;
    ratio = std::pow(ratio, 1./3.);
    if (ratio > 1)  { ratio = 1; }
    reductionRatio.set(ratio,ratio,ratio);
  }
}

//______________________________________________________________________________
void G4Voxelizer::BuildReduceVoxels(std::vector<G4double> boundaries[],
                                    G4ThreeVector reductionRatio)
{
  for (G4int k = 0; k <= 2; ++k)
  {
    std::vector<G4int> &candidatesCount = fCandidatesCounts[k];
    G4int max = candidatesCount.size();
    std::vector<G4VoxelInfo> voxels(max);
    G4VoxelComparator comp(voxels);
    std::set<G4int, G4VoxelComparator> voxelSet(comp);
    std::vector<G4int> mergings;

    for (G4int j = 0; j < max; ++j)
    {
      G4VoxelInfo &voxel = voxels[j];
      voxel.count = candidatesCount[j];
      voxel.previous = j - 1;
      voxel.next = j + 1;
      voxels[j] = voxel;
    }

    for (G4int j = 0; j < max - 1; ++j)  { voxelSet.insert(j); }
      // we go to size-1 to make sure we will not merge the last element

    G4double reduction = reductionRatio[k];
    if (reduction != 0)
    {
      G4int count = 0, currentCount;
      while ((currentCount = voxelSet.size()) > 2) 
      {
        G4double currentRatio = 1 - (G4double) count / max;
        if ((currentRatio <= reduction) && (currentCount <= 1000))
          break;
        const G4int pos = *voxelSet.begin();
        mergings.push_back(pos + 1);

        G4VoxelInfo& voxel = voxels[pos];
        G4VoxelInfo& nextVoxel = voxels[voxel.next];

        if (voxelSet.erase(pos) != 1)
        {
          ;// k = k;
        }
        if (voxel.next != max - 1)
          if (voxelSet.erase(voxel.next) != 1)
          {
            ;// k = k;
          }
        if (voxel.previous != -1)
          if (voxelSet.erase(voxel.previous) != 1)
          {
            ;// k = k;
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
      }  // Loop checking, 13.08.2015, G.Cosmo
    }

    if (mergings.size())
    {
      std::sort(mergings.begin(), mergings.end());

      const std::vector<G4double>& boundary = boundaries[k];
      int mergingsSize = mergings.size();
      vector<G4double> reducedBoundary;
      G4int skip = mergings[0], i = 0;
      max = boundary.size();
      for (G4int j = 0; j < max; ++j)
      {
        if (j != skip)
        {
          reducedBoundary.push_back(boundary[j]);
        }
        else if (++i < mergingsSize)
        {
          skip = mergings[i];
        }
      }
      boundaries[k] = reducedBoundary;
    }
/*
    G4int count = 0;
    while (true)    // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double reduction = reductionRatio[k];
      if (reduction == 0)
        break;
      G4int currentCount = voxelSet.size();
      if (currentCount <= 2)
        break;
      G4double currentRatio = 1 - (G4double) count / max;
      if (currentRatio <= reduction && currentCount <= 1000)
        break;
      const G4int pos = *voxelSet.begin();
      mergings.push_back(pos);

      G4VoxelInfo &voxel = voxels[pos];
      G4VoxelInfo &nextVoxel = voxels[voxel.next];

      voxelSet.erase(pos);
      if (voxel.next != max - 1) { voxelSet.erase(voxel.next); }
      if (voxel.previous != -1)  { voxelSet.erase(voxel.previous); }

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

    if (mergings.size())
    {
      std::sort(mergings.begin(), mergings.end());

      std::vector<G4double> &boundary = boundaries[k];
      std::vector<G4double> reducedBoundary(boundary.size() - mergings.size());
      G4int skip = mergings[0] + 1, cur = 0, i = 0;
      max = boundary.size();
      for (G4int j = 0; j < max; ++j)
      {
        if (j != skip)
        {
          reducedBoundary[cur++] = boundary[j];
        }
        else
        {
          if (++i < (G4int)mergings.size())  { skip = mergings[i] + 1; }
        }
      }
      boundaries[k] = reducedBoundary;
    }
*/
  }
}

//______________________________________________________________________________
void G4Voxelizer::BuildReduceVoxels2(std::vector<G4double> boundaries[],
                                     G4ThreeVector reductionRatio)
{
  for (G4int k = 0; k <= 2; ++k)
  {
    std::vector<G4int> &candidatesCount = fCandidatesCounts[k];
    G4int max = candidatesCount.size();
    G4int total = 0;
    for (G4int i = 0; i < max; ++i) total += candidatesCount[i];

    G4double reduction = reductionRatio[k];
    if (reduction == 0)
      break;

    G4int destination = (G4int) (reduction * max) + 1;
    if (destination > 1000) destination = 1000;
    if (destination < 2) destination = 2;
    G4double average = ((G4double)total / max) / reduction;

    std::vector<G4int> mergings;

    std::vector<G4double> &boundary = boundaries[k];
    std::vector<G4double> reducedBoundary(destination);

    G4int sum = 0, cur = 0;
    for (G4int i = 0; i < max; ++i)
    {
      sum += candidatesCount[i];
      if (sum > average * (cur + 1) || i == 0)
      {
        G4double val = boundary[i];
        reducedBoundary[cur] = val;
        cur++;
        if (cur == destination)
          break;
      }
    }
    reducedBoundary[destination-1] = boundary[max];
    boundaries[k] = reducedBoundary;
  }
}

//______________________________________________________________________________
void G4Voxelizer::Voxelize(std::vector<G4VSolid*>& solids,
                           std::vector<G4Transform3D>& transforms)
{
  BuildVoxelLimits(solids, transforms);
  BuildBoundaries();
  BuildBitmasks(fBoundaries, fBitmasks);
  BuildBoundingBox();
  BuildEmpty(); // this does not work well for multi-union,
                // actually only makes performance slower,
                // these are only pre-calculated but not used by multi-union

  for (G4int i = 0; i < 3; ++i)
  {
    fCandidatesCounts[i].resize(0);
  }
}

//______________________________________________________________________________
void G4Voxelizer::CreateMiniVoxels(std::vector<G4double> boundaries[],
                                   G4SurfBits bitmasks[])
{
  std::vector<G4int> voxel(3), maxVoxels(3);
  for (G4int i = 0; i <= 2; ++i) maxVoxels[i] = boundaries[i].size();

  G4ThreeVector point;
  for (voxel[2] = 0; voxel[2] < maxVoxels[2] - 1; ++voxel[2])
  {
    for (voxel[1] = 0; voxel[1] < maxVoxels[1] - 1; ++voxel[1])
    {
      for (voxel[0] = 0; voxel[0] < maxVoxels[0] - 1; ++voxel[0])
      {
        std::vector<G4int> candidates;
        if (GetCandidatesVoxelArray(voxel, bitmasks, candidates, 0))
        {
          // find a box for corresponding non-empty voxel
          G4VoxelBox box;
          for (G4int i = 0; i <= 2; ++i)
          {
            G4int index = voxel[i];
            const std::vector<G4double> &boundary = boundaries[i];
            G4double hlen = 0.5 * (boundary[index+1] - boundary[index]);
            box.hlen[i] = hlen;
            box.pos[i] = boundary[index] + hlen;
          }
          fVoxelBoxes.push_back(box);
          std::vector<G4int>(candidates).swap(candidates);
          fVoxelBoxesCandidates.push_back(candidates);
        }
      }
    }
  }
}

//______________________________________________________________________________
void G4Voxelizer::Voxelize(std::vector<G4VFacet*>& facets)
{
  G4int maxVoxels = fMaxVoxels;
  G4ThreeVector reductionRatio = fReductionRatio;

  G4int size = facets.size();
  if (size < 10)
  {
    for (G4int i = 0; i < (G4int) facets.size(); ++i)
    { 
      if (facets[i]->GetNumberOfVertices() > 3) size++;
    }
  }

  if ((size >= 10 || maxVoxels > 0) && maxVoxels != 0 && maxVoxels != 1)
  {
#ifdef G4SPECSDEBUG
    G4cout << "Building voxel limits..." << G4endl;
#endif
    
    BuildVoxelLimits(facets);

#ifdef G4SPECSDEBUG
    G4cout << "Building boundaries..." << G4endl;
#endif
    
    BuildBoundaries();

#ifdef G4SPECSDEBUG
    G4cout << "Building bitmasks..." << G4endl;
#endif
    
    BuildBitmasks(fBoundaries, 0, true);

    if (maxVoxels < 0 && reductionRatio == G4ThreeVector())
    {
      maxVoxels = fTotalCandidates;
      if (fTotalCandidates > 1000000) maxVoxels = 1000000;
    }

    SetReductionRatio(maxVoxels, reductionRatio);

    fCountOfVoxels = CountVoxels(fBoundaries);

#ifdef G4SPECSDEBUG
    G4cout << "Total number of voxels: " << fCountOfVoxels << G4endl;
#endif

    BuildReduceVoxels2(fBoundaries, reductionRatio);

    fCountOfVoxels = CountVoxels(fBoundaries);

#ifdef G4SPECSDEBUG
    G4cout << "Total number of voxels after reduction: "
           << fCountOfVoxels << G4endl;
#endif

#ifdef G4SPECSDEBUG
    G4cout << "Building bitmasks..." << G4endl;
#endif
    
    BuildBitmasks(fBoundaries, fBitmasks);

    G4ThreeVector reductionRatioMini;

    G4SurfBits bitmasksMini[3];

    // section for building mini voxels

    std::vector<G4double> miniBoundaries[3];

    for (G4int i = 0; i <= 2; ++i)  { miniBoundaries[i] = fBoundaries[i]; }

    G4int voxelsCountMini = (fCountOfVoxels >= 1000)
                          ? 100 : fCountOfVoxels / 10;

    SetReductionRatio(voxelsCountMini, reductionRatioMini);

#ifdef G4SPECSDEBUG
    G4cout << "Building reduced voxels..." << G4endl;
#endif
    
    BuildReduceVoxels(miniBoundaries, reductionRatioMini);

#ifdef G4SPECSDEBUG
    G4int total = CountVoxels(miniBoundaries);
    G4cout << "Total number of mini voxels: " << total << G4endl;
#endif

#ifdef G4SPECSDEBUG
    G4cout << "Building mini bitmasks..." << G4endl;
#endif
    
    BuildBitmasks(miniBoundaries, bitmasksMini);

#ifdef G4SPECSDEBUG
    G4cout << "Creating Mini Voxels..." << G4endl;
#endif
    
    CreateMiniVoxels(miniBoundaries, bitmasksMini);

#ifdef G4SPECSDEBUG
    G4cout << "Building bounding box..." << G4endl;
#endif
    
    BuildBoundingBox();

#ifdef G4SPECSDEBUG
    G4cout << "Building empty..." << G4endl;
#endif
    
    BuildEmpty();

#ifdef G4SPECSDEBUG
    G4cout << "Deallocating unnecessary fields during runtime..." << G4endl;
#endif    
    // deallocate fields unnecessary during runtime
    //
    fBoxes.resize(0);
    for (G4int i = 0; i < 3; ++i)
    {
      fCandidatesCounts[i].resize(0);
      fBitmasks[i].Clear();
    }
  }
}


//______________________________________________________________________________
void G4Voxelizer::GetCandidatesVoxel(std::vector<G4int> &voxels)
{
  // "GetCandidates" should compute which solids are possibly contained in
  // the voxel defined by the three slices characterized by the passed indexes.

  G4cout << "   Candidates in voxel [" << voxels[0] << " ; " << voxels[1]
         << " ; " << voxels[2] << "]: ";
  std::vector<G4int> candidates;
  G4int count = GetCandidatesVoxelArray(voxels, candidates);
  G4cout << "[ ";
  for (G4int i = 0; i < count; ++i) G4cout << candidates[i];
  G4cout << "]  " << G4endl;
}

//______________________________________________________________________________
void G4Voxelizer::FindComponentsFastest(unsigned int mask,
                                         std::vector<G4int> &list, G4int i)
{
  for (G4int byte = 0; byte < (G4int) (sizeof(unsigned int)); byte++)
  {
    if (G4int maskByte = mask & 0xFF)
    {
      for (G4int bit = 0; bit < 8; bit++)
      {
        if (maskByte & 1) 
          { list.push_back(8*(sizeof(unsigned int)*i+ byte) + bit); }
        if (!(maskByte >>= 1)) break;
      }
    }
    mask >>= 8;
  }
}

//______________________________________________________________________________
void G4Voxelizer::TransformLimits(G4ThreeVector& min, G4ThreeVector& max,
                                  const G4Transform3D& transformation) const
{
  // The goal of this method is to convert the quantities min and max
  // (representing the bounding box of a given solid in its local frame)
  // to the main frame, using "transformation"

  G4ThreeVector vertices[8] =     // Detemination of the vertices thanks to
  {                               // the extension of each solid:
    G4ThreeVector(min.x(), min.y(), min.z()), // 1st vertice:
    G4ThreeVector(min.x(), max.y(), min.z()), // 2nd vertice:
    G4ThreeVector(max.x(), max.y(), min.z()),
    G4ThreeVector(max.x(), min.y(), min.z()),
    G4ThreeVector(min.x(), min.y(), max.z()),
    G4ThreeVector(min.x(), max.y(), max.z()),
    G4ThreeVector(max.x(), max.y(), max.z()),
    G4ThreeVector(max.x(), min.y(), max.z())
  };

  min.set(kInfinity,kInfinity,kInfinity);
  max.set(-kInfinity,-kInfinity,-kInfinity);

  // Loop on th vertices
  G4int limit = sizeof(vertices) / sizeof(G4ThreeVector);
  for (G4int i = 0 ; i < limit; ++i)
  {
    // From local frame to the global one:
    // Current positions on the three axis:
    G4ThreeVector current = GetGlobalPoint(transformation, vertices[i]);

    // If need be, replacement of the min & max values:
    if (current.x() > max.x()) max.setX(current.x());
    if (current.x() < min.x()) min.setX(current.x());

    if (current.y() > max.y()) max.setY(current.y());
    if (current.y() < min.y()) min.setY(current.y());

    if (current.z() > max.z()) max.setZ(current.z());
    if (current.z() < min.z()) min.setZ(current.z());
  }
}

//______________________________________________________________________________
G4int G4Voxelizer::GetCandidatesVoxelArray(const G4ThreeVector &point,
                          std::vector<G4int> &list, G4SurfBits *crossed) const
{
  // Method returning the candidates corresponding to the passed point

  list.clear();

  for (G4int i = 0; i <= 2; ++i)
  {
    if(point[i] < fBoundaries[i].front() || point[i] >= fBoundaries[i].back()) 
      return 0;
  }

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
      G4int slice;
      if (fBoundaries[0].size() > 2)
      {
        slice = BinarySearch(fBoundaries[0], point.x());
        if (!(mask = ((unsigned int*) fBitmasks[0].fAllBits)[slice]))
          return 0;
      }
      if (fBoundaries[1].size() > 2)
      {
        slice = BinarySearch(fBoundaries[1], point.y());
        if (!(mask &= ((unsigned int*) fBitmasks[1].fAllBits)[slice]))
          return 0;
      }
      if (fBoundaries[2].size() > 2)
      {
        slice = BinarySearch(fBoundaries[2], point.z());
        if (!(mask &= ((unsigned int*) fBitmasks[2].fAllBits)[slice]))
          return 0;
      }
      if (crossed && (!(mask &= ~((unsigned int*)crossed->fAllBits)[0])))
        return 0;

      FindComponentsFastest(mask, list, 0);
    }
    else
    {
      unsigned int* masks[3], mask; // masks for X,Y,Z axis
      for (G4int i = 0; i <= 2; ++i)
      {
        G4int slice = BinarySearch(fBoundaries[i], point[i]);
        masks[i] = ((unsigned int*) fBitmasks[i].fAllBits)
                 + slice * fNPerSlice;
      }
      unsigned int* maskCrossed = crossed
                                ? (unsigned int*)crossed->fAllBits : 0;

      for (G4int i = 0 ; i < fNPerSlice; ++i)
      {
        // Logic "and" of the masks along the 3 axes x, y, z:
        // removing "if (!" and ") continue" => slightly slower
        //
        if (!(mask = masks[0][i])) continue;
        if (!(mask &= masks[1][i])) continue;
        if (!(mask &= masks[2][i])) continue;
        if (maskCrossed && !(mask &= ~maskCrossed[i])) continue;

        FindComponentsFastest(mask, list, i);
      }
    }
/*
    if (fNPerSlice == 1)
    {
      unsigned int mask;
      G4int slice = BinarySearch(fBoundaries[0], point.x()); 
      if (!(mask = ((unsigned int *) fBitmasks[0].fAllBits)[slice]
      )) return 0;
      slice = BinarySearch(fBoundaries[1], point.y());
      if (!(mask &= ((unsigned int *) fBitmasks[1].fAllBits)[slice]
      )) return 0;
      slice = BinarySearch(fBoundaries[2], point.z());
      if (!(mask &= ((unsigned int *) fBitmasks[2].fAllBits)[slice]
      )) return 0;
      if (crossed && (!(mask &= ~((unsigned int *)crossed->fAllBits)[0])))
        return 0;

      FindComponentsFastest(mask, list, 0);
    }
    else
    {
      unsigned int *masks[3], mask; // masks for X,Y,Z axis
      for (G4int i = 0; i <= 2; ++i)
      {
        G4int slice = BinarySearch(fBoundaries[i], point[i]); 
        masks[i] = ((unsigned int *) fBitmasks[i].fAllBits) + slice*fNPerSlice;
      }
      unsigned int *maskCrossed =
        crossed ? (unsigned int *)crossed->fAllBits : 0;

      for (G4int i = 0 ; i < fNPerSlice; ++i)
      {
        // Logic "and" of the masks along the 3 axes x, y, z:
        // removing "if (!" and ") continue" => slightly slower
        //
        if (!(mask = masks[0][i])) continue;
        if (!(mask &= masks[1][i])) continue;
        if (!(mask &= masks[2][i])) continue;
        if (maskCrossed && !(mask &= ~maskCrossed[i])) continue;

        FindComponentsFastest(mask, list, i);
      }
    }
*/
  }
  return list.size();
}

//______________________________________________________________________________
G4int
G4Voxelizer::GetCandidatesVoxelArray(const std::vector<G4int> &voxels,
                                     const G4SurfBits bitmasks[],
                                     std::vector<G4int> &list,
                                     G4SurfBits *crossed) const
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
      if (!(mask = ((unsigned int *) bitmasks[0].fAllBits)[voxels[0]]))
        return 0;
      if (!(mask &= ((unsigned int *) bitmasks[1].fAllBits)[voxels[1]]))
        return 0;
      if (!(mask &= ((unsigned int *) bitmasks[2].fAllBits)[voxels[2]]))
        return 0;
      if (crossed && (!(mask &= ~((unsigned int *)crossed->fAllBits)[0])))
        return 0;

      FindComponentsFastest(mask, list, 0);
    }
    else
    {
      unsigned int *masks[3], mask; // masks for X,Y,Z axis
      for (G4int i = 0; i <= 2; ++i)
      {
        masks[i] = ((unsigned int *) bitmasks[i].fAllBits)
                 + voxels[i]*fNPerSlice;
      }
      unsigned int *maskCrossed = crossed
                                ? (unsigned int *)crossed->fAllBits : 0;

      for (G4int i = 0 ; i < fNPerSlice; ++i)
      {
        // Logic "and" of the masks along the 3 axes x, y, z:
        // removing "if (!" and ") continue" => slightly slower
        //
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
G4int
G4Voxelizer::GetCandidatesVoxelArray(const std::vector<G4int> &voxels,
                          std::vector<G4int> &list, G4SurfBits *crossed) const
{
  // Method returning the candidates corresponding to the passed point

  return GetCandidatesVoxelArray(voxels, fBitmasks, list, crossed);
}

//______________________________________________________________________________
G4bool G4Voxelizer::Contains(const G4ThreeVector &point) const
{
  for (G4int i = 0; i < 3; ++i)
  {
    if (point[i] < fBoundaries[i].front() || point[i] > fBoundaries[i].back())
      return false;
  }
  return true;
}

//______________________________________________________________________________
G4double
G4Voxelizer::DistanceToFirst(const G4ThreeVector &point,
                             const G4ThreeVector &direction) const
{
  G4ThreeVector pointShifted = point - fBoundingBoxCenter;
  G4double shift = fBoundingBox.DistanceToIn(pointShifted, direction);
  return shift;
}

//______________________________________________________________________________
G4double
G4Voxelizer::DistanceToBoundingBox(const G4ThreeVector &point) const
{
  G4ThreeVector pointShifted = point - fBoundingBoxCenter;
  G4double shift = MinDistanceToBox(pointShifted, fBoundingBoxSize);
  return shift;
}

//______________________________________________________________________________
G4double
G4Voxelizer::MinDistanceToBox (const G4ThreeVector &aPoint,
                               const G4ThreeVector &f)
{
  // Estimates the isotropic safety from a point outside the current solid to
  // any of its surfaces. The algorithm may be accurate or should provide a
  // fast underestimate.

  G4double safe, safx, safy, safz;
  safe = safx = -f.x() + std::abs(aPoint.x());
  safy = -f.y() + std::abs(aPoint.y());
  if ( safy > safe ) safe = safy;
  safz = -f.z() + std::abs(aPoint.z());
  if ( safz > safe ) safe = safz;
  if (safe < 0.0) return 0.0; // point is inside

  G4double safsq = 0.0;
  G4int count = 0;
  if ( safx > 0 ) { safsq += safx*safx; count++; }
  if ( safy > 0 ) { safsq += safy*safy; count++; }
  if ( safz > 0 ) { safsq += safz*safz; count++; }
  if (count == 1) return safe;
  return std::sqrt(safsq);
}

//______________________________________________________________________________
G4double
G4Voxelizer::DistanceToNext(const G4ThreeVector &point,
                            const G4ThreeVector &direction,
                                  std::vector<G4int> &curVoxel) const
{
  G4double shift = kInfinity;

  G4int cur = 0; // the smallest index, which would be than increased
  for (int i = 0; i <= 2; ++i)
  {
    // Looking for the next voxels on the considered direction X,Y,Z axis
    //
    const std::vector<G4double>& boundary = fBoundaries[i];
    G4int index = curVoxel[i];
    if (direction[i] >= 1e-10)
    {
      ++index;
    }
    else
    {
      if (direction[i] > -1e-10)
        continue;
    }
    G4double dif = boundary[index] - point[i];
    G4double distance = dif / direction[i];

    if (shift > distance)
    {
      shift = distance;
      cur = i;
    }
  }

  if (shift != kInfinity)
  {
    // updating current voxel using the index corresponding
    // to the closest voxel boundary on the ray

    if (direction[cur] > 0)
    {
      if (++curVoxel[cur] >= (G4int) fBoundaries[cur].size() - 1)
        shift = kInfinity;
    }
    else
    {
      if (--curVoxel[cur] < 0)
        shift = kInfinity;
    }
  }

/*
  for (G4int i = 0; i <= 2; ++i)
  {
    // Looking for the next voxels on the considered direction X,Y,Z axis
    //
    const std::vector<G4double> &boundary = fBoundaries[i];
    G4int cur = curVoxel[i];
    if(direction[i] >= 1e-10)
    {
        if (boundary[++cur] - point[i] < fTolerance) // make sure shift would
        if (++cur >= (G4int) boundary.size())      // be non-zero
          continue;
    }
    else 
    {
      if(direction[i] <= -1e-10) 
      {
        if (point[i] - boundary[cur] < fTolerance) // make sure shift would
          if (--cur < 0)                           // be non-zero
            continue;
      }
      else 
        continue;
    }
    G4double dif = boundary[cur] - point[i];
    G4double distance = dif / direction[i];

    if (shift > distance) 
      shift = distance;
  }
*/
  return shift;
}

//______________________________________________________________________________
G4bool
G4Voxelizer::UpdateCurrentVoxel(const G4ThreeVector &point,
                                const G4ThreeVector &direction,
                                std::vector<G4int> &curVoxel) const
{
  for (G4int i = 0; i <= 2; ++i)
  {
    G4int index = curVoxel[i];
    const std::vector<G4double> &boundary = fBoundaries[i];

    if (direction[i] > 0) 
    {
      if (point[i] >= boundary[++index]) 
        if (++curVoxel[i] >= (G4int) boundary.size() - 1)
          return false;
    }
    else
    {
      if (point[i] < boundary[index]) 
        if (--curVoxel[i] < 0) 
          return false;
    }
#ifdef G4SPECSDEBUG
    G4int indexOK = BinarySearch(boundary, point[i]);
    if (curVoxel[i] != indexOK)
      curVoxel[i] = indexOK; // put breakpoint here
#endif
  }
  return true;
}

//______________________________________________________________________________
void G4Voxelizer::SetMaxVoxels(G4int max)
{
  fMaxVoxels = max;
  fReductionRatio.set(0,0,0);
}

//______________________________________________________________________________
void G4Voxelizer::SetMaxVoxels(const G4ThreeVector &ratioOfReduction)
{
  fMaxVoxels = -1;
  fReductionRatio = ratioOfReduction;
}

//______________________________________________________________________________
void G4Voxelizer::SetDefaultVoxelsCount(G4int count)
{
  fDefaultVoxelsCount = count;
}

//______________________________________________________________________________
G4int G4Voxelizer::GetDefaultVoxelsCount()
{
  return fDefaultVoxelsCount;
}

//______________________________________________________________________________
G4int G4Voxelizer::AllocatedMemory()
{
  G4int size = fEmpty.GetNbytes();
  size += fBoxes.capacity() * sizeof(G4VoxelBox);
  size += sizeof(G4double) * (fBoundaries[0].capacity()
        + fBoundaries[1].capacity() + fBoundaries[2].capacity());
  size += sizeof(G4int) * (fCandidatesCounts[0].capacity()
        + fCandidatesCounts[1].capacity() + fCandidatesCounts[2].capacity());
  size += fBitmasks[0].GetNbytes() + fBitmasks[1].GetNbytes()
        + fBitmasks[2].GetNbytes();

  G4int csize = fCandidates.size();
  for (G4int i = 0; i < csize; ++i)
  {
    size += sizeof(vector<G4int>) + fCandidates[i].capacity() * sizeof(G4int);
  }

  return size;
}
