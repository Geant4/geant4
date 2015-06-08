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
// UTessellatedSolid
//
// 11.07.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UTessellatedSolid.hh"

#include <algorithm>
#include <stack>
#include <list>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

/// TODO: make a benchmark for automatic selection of number of voxels. random voxels will be selected,
/// than for them methods distancetoin/out and inside will be launched. eventually, find out from Geant4 how it is done

///////////////////////////////////////////////////////////////////////////////
//
// Standard contructor has blank name and defines no fFacets.
//

void UTessellatedSolid::Initialize()
{
  fgToleranceHalf = 0.5 * fgTolerance;

  // I recommend keeping NULL here instead of 0. c++11 provides nullptr, than it could be easily replaced everywhere using simple case sensitive, whole words replacement in all files

  fCubicVolume = 0.;
  fSurfaceArea = 0.;

  fGeometryType = "TessellatedSolid";
  fSolidClosed  = false;

  fMinExtent.Set(UUtils::kInfinity);
  fMaxExtent.Set(-UUtils::kInfinity);

  SetRandomVectors();
}

UTessellatedSolid::UTessellatedSolid() : VUSolid("dummy")
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor. Simple define name and geometry type - no fFacets
// to detine.
//
UTessellatedSolid::UTessellatedSolid(const string& name)
  : VUSolid(name)
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UTessellatedSolid::UTessellatedSolid(__void__&) : VUSolid("")
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
UTessellatedSolid::~UTessellatedSolid()
{
  DeleteObjects();
}

///////////////////////////////////////////////////////////////////////////////
//
// Define copy constructor.
//
UTessellatedSolid::UTessellatedSolid(const UTessellatedSolid& ss)
  : VUSolid(ss)
{
  Initialize();

  CopyObjects(ss);
}

///////////////////////////////////////////////////////////////////////////////
//
// Define assignment operator.
//
UTessellatedSolid& UTessellatedSolid::operator= (const UTessellatedSolid& ss)
{
  if (&ss == this) return *this;

  // Copy base class data
  VUSolid::operator=(ss);

  DeleteObjects();

  Initialize();

  CopyObjects(ss);

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::DeleteObjects()
{
  int size = fFacets.size();
  for (int i = 0; i < size; ++i)
    delete fFacets[i];
  fFacets.clear();
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::CopyObjects(const UTessellatedSolid& ss)
{
  UVector3 reductionRatio;
  int fmaxVoxels = fVoxels.GetMaxVoxels(reductionRatio);
  if (fmaxVoxels < 0)
    fVoxels.SetMaxVoxels(reductionRatio);
  else
    fVoxels.SetMaxVoxels(fmaxVoxels);

  int n = ss.GetNumberOfFacets();
  for (int i = 0; i < n; ++i)
  {
    VUFacet* facetClone = (ss.GetFacet(i))->GetClone();
    AddFacet(facetClone);
  }
  if (ss.GetSolidClosed()) SetSolidClosed(true);
}


///////////////////////////////////////////////////////////////////////////////
//
// Add a facet to the facet list.  Note that you can add, but you cannot
// delete.
//
bool UTessellatedSolid::AddFacet(VUFacet* aFacet)
{
  // Add the facet to the vector.
  if (fSolidClosed)
  {
    UUtils::Exception("UTessellatedSolid::AddFacet()", "GeomSolids1002",
      UWarning, 1, "Attempt to add facets when solid is closed.");
    return false;
  }
  else if (aFacet->IsDefined())
  {
    set<UVertexInfo, UVertexComparator>::iterator begin = fFacetList.begin(), end = fFacetList.end(), pos, it;
    UVector3 p = aFacet->GetCircumcentre();
    UVertexInfo value;
    value.id = fFacetList.size();
    value.mag2 = p.x() + p.y() + p.z();

    bool found = false;
    if (!OutsideOfExtent(p, fgTolerance))
    {
      double fgTolerance3 = 3 * fgTolerance;
      pos = fFacetList.lower_bound(value);

      it = pos;
      while (!found && it != end)
      {
        int id = (*it).id;
        VUFacet* facet = fFacets[id];
        UVector3 q = facet->GetCircumcentre();
        found = (facet == aFacet);
        if (found) break;
        double dif = q.x() + q.y() + q.z() - value.mag2;
        if (dif > fgTolerance3) break;
        it++;
      }

      if (fFacets.size() > 1)
      {
        it = pos;
        while (!found && it != begin)
        {
          --it;
          int id = (*it).id;
          VUFacet* facet = fFacets[id];
          UVector3 q = facet->GetCircumcentre();
          found = (facet == aFacet);
          if (found) break;
          double dif = q.x() + q.y() + q.z() - q.Mag2();
          if (dif > fgTolerance3) break;
        }
      }
    }

    if (!found)
    {
      fFacets.push_back(aFacet);
      fFacetList.insert(value);
    }
    return true;
  }
  else
  {
    UUtils::Exception("UTessellatedSolid::AddFacet()", "GeomSolids1002",
      UWarning, 1, "Attempt to add facet not properly defined.");
    aFacet->StreamInfo(cout);
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
int UTessellatedSolid::SetAllUsingStack(const vector<int>& voxel, const vector<int>& max, bool status, UBits& checked)
{
  vector<int> xyz = voxel;
  stack<vector<int> > pos;
  pos.push(xyz);
  int filled = 0;
  int cc = 0, nz = 0;

  vector<int> candidates;

  while (!pos.empty())
  {
    xyz = pos.top();
    pos.pop();
    int index = fVoxels.GetVoxelsIndex(xyz);
    if (!checked[index])
    {
      checked.SetBitNumber(index, true);
      cc++;

      if (fVoxels.IsEmpty(index))
      {
        filled++;

        fInsides.SetBitNumber(index, status);

        for (int i = 0; i <= 2; ++i)
        {
          if (xyz[i] < max[i] - 1)
          {
            xyz[i]++;
            pos.push(xyz);
            xyz[i]--;
          }

          if (xyz[i] > 0)
          {
            xyz[i]--;
            pos.push(xyz);
            xyz[i]++;
          }
        }
      }
      else
      {
        nz++;
      }
    }
  }
  return filled;
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::PrecalculateInsides()
{
  vector<int> voxel(3), maxVoxels(3);
  for (int i = 0; i <= 2; ++i) maxVoxels[i] = fVoxels.GetBoundary(i).size();
  int size = maxVoxels[0] * maxVoxels[1] * maxVoxels[2];

  UBits checked(size - 1);
  fInsides.Clear();
  fInsides.ResetBitNumber(size - 1);

  UVector3 point;
  for (voxel[2] = 0; voxel[2] < maxVoxels[2] - 1; ++voxel[2])
  {
    for (voxel[1] = 0; voxel[1] < maxVoxels[1] - 1; ++voxel[1])
    {
      for (voxel[0] = 0; voxel[0] < maxVoxels[0] - 1; ++voxel[0])
      {
        int index = fVoxels.GetVoxelsIndex(voxel);
        if (!checked[index] && fVoxels.IsEmpty(index))
        {
          for (int i = 0; i <= 2; ++i) point[i] = fVoxels.GetBoundary(i)[voxel[i]];
          bool inside = (bool)(InsideNoVoxels(point) == eInside);
          SetAllUsingStack(voxel, maxVoxels, inside, checked);
        }
        else checked.SetBitNumber(index);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::Voxelize()
{
#ifdef USPECSDEBUG
  cout << "Voxelizing...\n";
#endif
  fVoxels.Voxelize(fFacets);

  if (fVoxels.Empty().GetNbits())
  {
#ifdef USPECSDEBUG
    cout << "Precalculating Insides...\n";
#endif
    PrecalculateInsides();
  }
  fVoxels.BuildBoundingBox(fMinExtent, fMaxExtent, fgTolerance);
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::SetExtremeFacets()
{
  //
  // Compute extremeFacets, i.e. find those facets that have surface
  // planes that bound the volume.
  // Note that this is going to reject concaved surfaces as being extreme.  Also
  // note that if the vertex is on the facet, displacement is zero, so IsInside
  // returns true.  So will this work??  Need non-equality
  // "bool inside = displacement < 0.0;"
  // or
  // "bool inside = displacement <= -0.5*fgTolerance"
  // (Notes from PT 13/08/2007).
  //
  int size = fFacets.size();
  int vsize = fVertexList.size();

  for (int j = 0; j < size; ++j)
  {
    VUFacet& facet = *fFacets[j];

    bool isExtreme = true;
    for (int i = 0; i < vsize; ++i)
    {
      if (!facet.IsInside(fVertexList[i]))
      {
        isExtreme = false;
        break;
      }
    }
    if (isExtreme) fExtremeFacets.insert(&facet);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::CreateVertexList()
{
  // the algorithm will be:
  // we will have additional vertexListSorted, where all the items will be sorted by magnitude of vertice vector
  // new candidate for fVertexList - we will determine the position fo first item which would be within it'ss magnitude - 0.5*fgTolerance. we will go trough until we will reach > +0.5 fgTolerance. comparison (q-p).Mag() < 0.5*fgTolerance will be made
  // they can be just stored in std::vector, with custom insertion based on binary search

  set<UVertexInfo, UVertexComparator> vertexListSorted;
  set<UVertexInfo, UVertexComparator>::iterator begin = vertexListSorted.begin(), end = vertexListSorted.end(), pos, it;
  UVector3 p;
  UVertexInfo value;

  fVertexList.clear();
  int size = fFacets.size();
  vector<int> newIndex(100);

  double fgTolerance24 = fgTolerance * fgTolerance / 4.0;
  double fgTolerance3 = 3 * fgTolerance;
  for (int k = 0; k < size; ++k)
  {
    VUFacet& facet = *fFacets[k];
    int max = facet.GetNumberOfVertices();

    for (int i = 0; i < max; ++i)
    {
      p = facet.GetVertex(i);
      value.id = fVertexList.size();
      value.mag2 = p.x() + p.y() + p.z();

      bool found = false;
      int id = 0;
      if (!OutsideOfExtent(p, fgTolerance))
      {
        pos = vertexListSorted.lower_bound(value);

        it = pos;
        while (it != end)
        {
          id = (*it).id;
          UVector3 q = fVertexList[id];
          double dif = (q - p).Mag2();
          found = (dif < fgTolerance24);
          if (found) break;
          dif = q.x() + q.y() + q.z() - value.mag2;
          if (dif > fgTolerance3) break;
          it++;
        }

        if (!found && fVertexList.size() > 1)
        {
          it = pos;
          while (it != begin)
          {
            --it;
            id = (*it).id;
            UVector3 q = fVertexList[id];
            double dif = (q - p).Mag2();
            found = (dif < fgTolerance24);
            if (found) break;
            dif = value.mag2 - (q.x() + q.y() + q.z());
            if (dif > fgTolerance) break;
          }
        }
      }

      //    cout << "Total checked: " << checked << " from " << fVertexList.size() << endl;

      if (!found)
      {
#ifdef USPECSDEBUG
        cout << p.x() << ":" << p.y() << ":" << p.z() << endl;
        cout << "Adding new vertex #" << i << " of facet " << k << " id " << value.id << endl;
        cout << "===" << endl;
#endif
        fVertexList.push_back(p);
        vertexListSorted.insert(value);
        begin = vertexListSorted.begin();
        end = vertexListSorted.end();
        newIndex[i] = value.id;

        //
        // Now update the maximum x, y and z limits of the volume.
        //
        if (value.id == 0) fMinExtent = fMaxExtent = p;
        else
        {
          if (p.x() > fMaxExtent.x()) fMaxExtent.x() = p.x();
          else if (p.x() < fMinExtent.x()) fMinExtent.x() = p.x();
          if (p.y() > fMaxExtent.y()) fMaxExtent.y() = p.y();
          else if (p.y() < fMinExtent.y()) fMinExtent.y() = p.y();
          if (p.z() > fMaxExtent.z()) fMaxExtent.z() = p.z();
          else if (p.z() < fMinExtent.z()) fMinExtent.z() = p.z();
        }
      }
      else
      {
#ifdef USPECSDEBUG
        cout << p.x() << ":" << p.y() << ":" << p.z() << endl;
        cout << "Vertex #" << i << " of facet " << k << " found, redirecting to " << id << endl;
        cout << "===" << endl;
#endif
        newIndex[i] = id;
      }
    }
    // only now it is possible to change vertices pointer
    facet.SetVertices(&fVertexList);
    for (int i = 0; i < max; i++)
      facet.SetVertexIndex(i, newIndex[i]);

  }
  vector<UVector3>(fVertexList).swap(fVertexList); // fVertexList.shrink_to_fit();

#ifdef DEBUG
  double previousValue = 0;
  for (res = vertexListSorted.begin(); res != vertexListSorted.end(); ++res)
  {
    int id = (*res).id;
    UVector3 vec = fVertexList[id];
    double value = abs(vec.Mag());
    if (previousValue > value)
      cout << "Error!" << "\n";
    previousValue = value;
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::DisplayAllocatedMemory()
{
  int without = AllocatedMemoryWithoutVoxels();
  int with = AllocatedMemory();
  double ratio = (double) with / without;
  cout << "Allocated memory without voxel overhead " << without << "; with " << with << "; ratio: " << ratio << endl;
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::SetSolidClosed(const bool t)
{
  if (t)
  {
    CreateVertexList();

    SetExtremeFacets();

    Voxelize();

#ifdef USPECSDEBUG
    DisplayAllocatedMemory();
#endif

  }
  fSolidClosed = t;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetSolidClosed
//
// Used to determine whether the solid is closed to adding further fFacets.
//
bool UTessellatedSolid::GetSolidClosed() const
{
  return fSolidClosed;
}

///////////////////////////////////////////////////////////////////////////////
//
// operator+=
//
// This operator allows the user to add two tessellated solids together, so
// that the solid on the left then includes all of the facets in the solid
// on the right.  Note that copies of the facets are generated, rather than
// using the original facet set of the solid on the right.
//
UTessellatedSolid& UTessellatedSolid::operator+=
(const UTessellatedSolid& right)
{
  int size = right.GetNumberOfFacets();
  for (int i = 0; i < size; ++i)
    AddFacet(right.GetFacet(i)->GetClone());

  return *this;
}


///////////////////////////////////////////////////////////////////////////////
//
// GetNumberOfFacets
//
int UTessellatedSolid::GetNumberOfFacets() const
{
  return fFacets.size();
}

///////////////////////////////////////////////////////////////////////////////
//
VUSolid::EnumInside UTessellatedSolid::InsideVoxels(const UVector3& p) const
{
  //
  // First the simple test - check if we're outside of the X-Y-Z extremes
  // of the tessellated solid.
  //
  if (OutsideOfExtent(p, fgTolerance))
    return eOutside;

  vector<int> startingVoxel(3);
  fVoxels.GetVoxel(startingVoxel, p);

  const vector<int>& startingCandidates = fVoxels.GetCandidates(startingVoxel);
  int limit = startingCandidates.size();
//  int limit = fVoxels.GetCandidatesVoxelArray(p, candidates, NULL);
  if (limit == 0 && fInsides.GetNbits())
  {
//    int index = fVoxels.GetPointIndex(p);
    int index = fVoxels.GetVoxelsIndex(startingVoxel);
    EnumInside location = fInsides[index] ? eInside : eOutside;
    return location;
  }

  double minDist = UUtils::kInfinity;

  for (int i = 0; i < limit; ++i)
  {
    int candidate = startingCandidates[i];
    VUFacet& facet = *fFacets[candidate];
    double dist = facet.Distance(p, minDist);
    if (dist < minDist) minDist = dist;
    if (dist <= fgToleranceHalf)
      return eSurface;
  }
  //
  //
  // The following is something of an adaptation of the method implemented by
  // Rickard Holmberg augmented with information from Schneider & Eberly,
  // "Geometric Tools for Computer Graphics," pp700-701, 2003.  In essence, we're
  // trying to determine whether we're inside the volume by projecting a few rays
  // and determining if the first surface crossed is has a normal vector between
  // 0 to pi/2 (out-going) or pi/2 to pi (in-going).  We should also avoid rays
  // which are nearly within the plane of the tessellated surface, and therefore
  // produce rays randomly.  For the moment, this is a bit over-engineered
  // (belt-braces-and-ducttape).
  //
  double distOut          = UUtils::kInfinity;
  double distIn           = UUtils::kInfinity;
  double distO            = 0.0;
  double distI            = 0.0;
  double distFromSurfaceO = 0.0;
  double distFromSurfaceI = 0.0;
  UVector3 normalO, normalI;
  bool crossingO          = false;
  bool crossingI          = false;
  VUSolid::EnumInside location          = VUSolid::eOutside;
//  VUSolid::EnumInside locationprime     = VUSolid::eOutside;
  int sm                   = 0;
  double shift;

  //  for (int i=0; i<nTry; ++i)
  //  {
  bool nearParallel = false;
  do
  {
    //
    //
    // We loop until we find direction where the vector is not nearly parallel
    // to the surface of any facet since this causes ambiguities.  The usual
    // case is that the angles should be sufficiently different, but there are 20
    // random directions to select from - hopefully sufficient.
    //
    distOut = distIn = UUtils::kInfinity;
    const UVector3& v = fRandir[sm];
    sm++;

    // This code could be voxelized by same algorithm, which is used for DistanceToOut.
    // we will traverse through fVoxels. we will call intersect only for those, which would be candidates
    // and was not checked before.

    //  double minDistance = UUtils::kInfinity;
//    UVector3 currentPoint = p;
    UVector3 direction = v.Unit();
//    UBits exclusion(fVoxels.GetBitsPerSlice());
    vector<int> curVoxel(3);
    curVoxel = startingVoxel;
//    double shiftBonus = VUSolid::Tolerance();

    bool crossed = false;
    bool started = true;
//    set<int> already;

    do
    {
      const vector<int>& candidates = started ? startingCandidates : fVoxels.GetCandidates(curVoxel);
      started = false;
      if (int candidatesCount = candidates.size())
      {
//          int candidatesCount = candidates.size();
        for (int i = 0 ; i < candidatesCount; ++i)
        {
          //            bits.SetBitNumber(candidate);
          int candidate = candidates[i];
          VUFacet& facet = *fFacets[candidate];

          crossingO = facet.Intersect(p, v, true, distO, distFromSurfaceO, normalO);
          crossingI = facet.Intersect(p, v, false, distI, distFromSurfaceI, normalI);

          if (crossingO || crossingI)
          {
            crossed = true;

//            double dot = std::fabs(normalO.Dot(v));
            nearParallel = (crossingO && std::fabs(normalO.Dot(v)) < dirTolerance) ||
                           (crossingI && std::fabs(normalI.Dot(v)) < dirTolerance);
            if (!nearParallel)
            {
              if (crossingO && distO > 0.0 && distO < distOut)
                distOut = distO;
              if (crossingI && distI > 0.0 && distI < distIn)
                distIn  = distI;
            }
            else break;
          }
        }
        if (nearParallel) break;
      }
      else
      {
        if (!crossed)
        {
          int index = fVoxels.GetVoxelsIndex(curVoxel);
          bool inside = fInsides[index];
          location = inside ? eInside : eOutside;
          // cout << (inside ? "I" : "O");
          return location;
        }
      }
      shift = fVoxels.DistanceToNext(p, direction, curVoxel);
    }
    while (shift != UUtils::kInfinity);

  }
  while (nearParallel && sm != fMaxTries);

  //
  //
  // Here we loop through the facets to find out if there is an intersection
  // between the ray and that facet.  The test if performed separately whether
  // the ray is entering the facet or exiting.
  //

#ifdef UVERBOSE
  if (sm == fMaxTries)
  {
    //
    // We've run out of random vector directions.  If nTries is set sufficiently
    // low (nTries <= 0.5*maxTries) then this would indicate that there is
    // something wrong with geometry.
    //
    std::ostringstream message;
    int oldprc = message.precision(16);
    message << "Cannot determine whether point is inside or outside volume!"
            << endl
            << "Solid name       = " << GetName()  << endl
            << "Geometry Type    = " << fGeometryType  << endl
            << "Number of facets = " << fFacets.size() << endl
            << "Position:"  << endl << endl
            << "p.x() = "   << p.x() / mm << " mm" << endl
            << "p.y() = "   << p.y() / mm << " mm" << endl
            << "p.z() = "   << p.z() / mm << " mm";
    message.precision(oldprc);
    UUtils::Exception("UTessellatedSolid::Inside()",
                      "GeomSolids1002", UWarning, 1, message.str().c_str());
  }
#endif
  //
  //
  // In the next if-then-elseif string the logic is as follows:
  // (1) You don't hit anything so cannot be inside volume, provided volume
  //     constructed correctly!
  // (2) Distance to inside (ie. nearest facet such that you enter facet) is
  //     shorter than distance to outside (nearest facet such that you exit
  //     facet) - on condition of safety distance - therefore we're outside.
  // (3) Distance to outside is shorter than distance to inside therefore we're
  //     inside.
  //
  if (distIn == UUtils::kInfinity && distOut == UUtils::kInfinity)
    location = eOutside;
  else if (distIn <= distOut - fgToleranceHalf)
    location = eOutside;
  else if (distOut <= distIn - fgToleranceHalf)
    location = eInside;
  // }

  // cout << " => " << (location == eInside ? "I" : "O") << endl;

  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
// VUSolid::EnumInside UTessellatedSolid::Inside (const UVector3 &p) const
//
// This method must return:
//    * kOutside if the point at offset p is outside the shape
//      boundaries plus fgTolerance/2,
//    * kSurface if the point is <= fgTolerance/2 from a surface, or
//    * kInside otherwise.
//
VUSolid::EnumInside UTessellatedSolid::Inside(const UVector3& aPoint) const
{
  VUSolid::EnumInside location;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    location = InsideVoxels(aPoint);
  }
  else
  {
    location = InsideNoVoxels(aPoint);
  }
#ifdef DEBUG
  if (location != insideNoVoxels)
    location = insideNoVoxels; // you can place a breakpoint here
#endif
  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
VUSolid::EnumInside UTessellatedSolid::InsideNoVoxels(const UVector3& p) const
{
  //
  // First the simple test - check if we're outside of the X-Y-Z extremes
  // of the tessellated solid.
  //
  if (OutsideOfExtent(p, fgTolerance))
    return eOutside;

  double minDist = UUtils::kInfinity;
  //
  //
  // Check if we are close to a surface
  //
  int size = fFacets.size();
  for (int i = 0; i < size; ++i)
  {
    VUFacet& facet = *fFacets[i];
    double dist = facet.Distance(p, minDist);
    if (dist < minDist) minDist = dist;
    if (dist <= fgToleranceHalf)
    {
      return eSurface;
    }
  }
  //
  //
  // The following is something of an adaptation of the method implemented by
  // Rickard Holmberg augmented with information from Schneider & Eberly,
  // "Geometric Tools for Computer Graphics," pp700-701, 2003.  In essence, we're
  // trying to determine whether we're inside the volume by projecting a few rays
  // and determining if the first surface crossed is has a normal vector between
  // 0 to pi/2 (out-going) or pi/2 to pi (in-going).  We should also avoid rays
  // which are nearly within the plane of the tessellated surface, and therefore
  // produce rays randomly.  For the moment, this is a bit over-engineered
  // (belt-braces-and-ducttape).
  //
#if USPECSDEBUG
  int nTry                = 7;
#else
  int nTry                = 3;
#endif
  double distOut          = UUtils::kInfinity;
  double distIn           = UUtils::kInfinity;
  double distO            = 0.0;
  double distI            = 0.0;
  double distFromSurfaceO = 0.0;
  double distFromSurfaceI = 0.0;
  UVector3 normalO(0.0, 0.0, 0.0);
  UVector3 normalI(0.0, 0.0, 0.0);
  bool crossingO          = false;
  bool crossingI          = false;
  VUSolid::EnumInside location          = VUSolid::eOutside;
  VUSolid::EnumInside locationprime     = VUSolid::eOutside;
  int sm = 0;

  for (int i = 0; i < nTry; ++i)
  {
    bool nearParallel = false;
    do
    {
      //
      //
      // We loop until we find direction where the vector is not nearly parallel
      // to the surface of any facet since this causes ambiguities.  The usual
      // case is that the angles should be sufficiently different, but there are 20
      // random directions to select from - hopefully sufficient.
      //
      distOut =  distIn = UUtils::kInfinity;
      UVector3 v = fRandir[sm];
      sm++;
      vector<VUFacet*>::const_iterator f = fFacets.begin();

      do
      {
        //
        //
        // Here we loop through the facets to find out if there is an intersection
        // between the ray and that facet.  The test if performed separately whether
        // the ray is entering the facet or exiting.
        //
        crossingO = ((*f)->Intersect(p, v, true, distO, distFromSurfaceO, normalO));
        crossingI = ((*f)->Intersect(p, v, false, distI, distFromSurfaceI, normalI));
        if (crossingO || crossingI)
        {
          nearParallel = (crossingO && std::fabs(normalO.Dot(v)) < dirTolerance) ||
                         (crossingI && std::fabs(normalI.Dot(v)) < dirTolerance);
          if (!nearParallel)
          {
            if (crossingO && distO > 0.0 && distO < distOut) distOut = distO;
            if (crossingI && distI > 0.0 && distI < distIn)  distIn  = distI;
          }
        }
      }
      while (!nearParallel && ++f != fFacets.end());
    }
    while (nearParallel && sm != fMaxTries);

#ifdef UVERBOSE
    if (sm == fMaxTries)
    {
      //
      //
      // We've run out of random vector directions.  If nTries is set sufficiently
      // low (nTries <= 0.5*maxTries) then this would indicate that there is
      // something wrong with geometry.
      //
      std::ostringstream message;
      int oldprc = message.precision(16);
      message << "Cannot determine whether point is inside or outside volume!"
              << endl
              << "Solid name       = " << GetName()  << endl
              << "Geometry Type    = " << fGeometryType  << endl
              << "Number of facets = " << fFacets.size() << endl
              << "Position:"  << endl << endl
              << "p.x() = "   << p.x() / mm << " mm" << endl
              << "p.y() = "   << p.y() / mm << " mm" << endl
              << "p.z() = "   << p.z() / mm << " mm";
      message.precision(oldprc);
      UUtils::Exception("UTessellatedSolid::Inside()",
                        "GeomSolids1002", UWarning, 1, message.str().c_str());
    }
#endif
    //
    //
    // In the next if-then-elseif string the logic is as follows:
    // (1) You don't hit anything so cannot be inside volume, provided volume
    //     constructed correctly!
    // (2) Distance to inside (ie. nearest facet such that you enter facet) is
    //     shorter than distance to outside (nearest facet such that you exit
    //     facet) - on condition of safety distance - therefore we're outside.
    // (3) Distance to outside is shorter than distance to inside therefore we're
    //     inside.
    //
    if (distIn == UUtils::kInfinity && distOut == UUtils::kInfinity)
      locationprime = eOutside;
    else if (distIn <= distOut - fgToleranceHalf)
      locationprime = eOutside;
    else if (distOut <= distIn - fgToleranceHalf)
      locationprime = eInside;

    if (i == 0) location = locationprime;
  }

  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.
//
bool UTessellatedSolid::Normal(const UVector3& p, UVector3& aNormal) const
{
  double minDist;
  VUFacet* facet = NULL;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    vector<int> curVoxel(3);
    fVoxels.GetVoxel(curVoxel, p);
    const vector<int>& candidates = fVoxels.GetCandidates(curVoxel);
//      fVoxels.GetCandidatesVoxelArray(p, candidates, NULL);

    if (int limit = candidates.size())
    {
      minDist = UUtils::kInfinity;
      for (int i = 0 ; i < limit ; ++i)
      {
        int candidate = candidates[i];
        VUFacet& f = *fFacets[candidate];
        double dist = f.Distance(p, minDist);
        if (dist < minDist) minDist = dist;
        if (dist <= fgToleranceHalf)
        {
          aNormal = f.GetSurfaceNormal();
          return true;
        }
      }
    }
    minDist = MinDistanceFacet(p, true, facet);
  }
  else
  {
    minDist = UUtils::kInfinity;
    int size = fFacets.size();
    for (int i = 0; i < size; ++i)
    {
      VUFacet& f = *fFacets[i];
      double dist = f.Distance(p, minDist);
      if (dist < minDist)
      {
        minDist  = dist;
        facet = &f;
      }
    }
  }

  if (minDist != UUtils::kInfinity)
  {
    if (facet) aNormal = facet->GetSurfaceNormal();
    return minDist <= fgToleranceHalf;
  }
  else
  {
#ifdef UVERBOSE
    std::ostringstream message;
    message << "Point p is not on surface !?" << endl
            << "          No facets found for point: " << p << " !" << endl
            << "          Returning approximated value for normal.";

    UUtils::Exception("UTessellatedSolid::SurfaceNormal(p)", "GeomSolids1002",
                      UWarning, 1, message.str().c_str());
#endif
    aNormal = (p.z() > 0 ? UVector3(0, 0, 1) : UVector3(0, 0, -1));
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// double DistanceToIn(const UVector3& p, const UVector3& v)
//
// Return the distance along the normalised vector v to the shape,
// from the point at offset p. If there is no intersection, return
// UUtils::kInfinity. The first intersection resulting from ‘leaving’ a
// surface/volume is discarded. Hence, this is tolerant of points on
// surface of shape.
//
double UTessellatedSolid::DistanceToInNoVoxels(const UVector3& p,
                                               const UVector3& v, double /*aPstep*/) const
{
  double minDist         = UUtils::kInfinity;
  double dist            = 0.0;
  double distFromSurface = 0.0;
  UVector3 normal;

#if USPECSDEBUG
  if (Inside(p) == kInside)
  {
    std::ostringstream message;
    int oldprc = message.precision(16) ;
    message << "Point p is already inside!?" << endl
            << "Position:"  << endl << endl
            << "   p.x() = "   << p.x() / mm << " mm" << endl
            << "   p.y() = "   << p.y() / mm << " mm" << endl
            << "   p.z() = "   << p.z() / mm << " mm" << endl
            << "DistanceToOut(p) == " << DistanceToOut(p);
    message.precision(oldprc) ;
    UUtils::Exception("UTriangularFacet::DistanceToIn(p,v)", "GeomSolids1002",
                      UWarning, 1, message.str().c_str());
  }
#endif

  int size = fFacets.size();
  for (int i = 0; i < size; ++i)
  {
    VUFacet& facet = *fFacets[i];
    if (facet.Intersect(p, v, false, dist, distFromSurface, normal))
    {
      //
      // Set minDist to the new distance to current facet if distFromSurface is in
      // positive direction and point is not at surface.  If the point is within
      // 0.5*fgTolerance of the surface, then force distance to be zero and
      // leave member function immediately (for efficiency), as proposed by & credit
      // to Akira Okumura.
      //
      if (distFromSurface > fgToleranceHalf && dist >= 0.0 && dist < minDist)
      {
        minDist  = dist;
      }
      if (-fgToleranceHalf <= dist && dist <= fgToleranceHalf)
      {
        return 0.0;
      }
      else if (distFromSurface > - fgToleranceHalf && distFromSurface < fgToleranceHalf) minDist = dist;

    }

  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::DistanceToOutNoVoxels(const UVector3& p, const UVector3& v, UVector3& aNormalVector, bool& aConvex, double /*aPstep*/) const
{
  double minDist         = UUtils::kInfinity;
  double dist            = 0.0;
  double distFromSurface = 0.0;
  UVector3 normal, minNormal;

#if USPECSDEBUG
  if (Inside(p) == kOutside)
  {
    std::ostringstream message;
    int oldprc = message.precision(16) ;
    message << "Point p is already outside!?" << endl
            << "Position:"  << endl << endl
            << "   p.x() = "   << p.x() / mm << " mm" << endl
            << "   p.y() = "   << p.y() / mm << " mm" << endl
            << "   p.z() = "   << p.z() / mm << " mm" << endl
            << "DistanceToIn(p) == " << DistanceToIn(p);
    message.precision(oldprc) ;
    UUtils::Exception("UTriangularFacet::DistanceToOut(p)", "GeomSolids1002",
                      UWarning, 1, message.str().c_str());
  }
#endif

  bool isExtreme = false;
  int size = fFacets.size();
  for (int i = 0; i < size; ++i)
  {
    VUFacet& facet = *fFacets[i];
    if (facet.Intersect(p, v, true, dist, distFromSurface, normal))
    {
      if (distFromSurface > 0.0 && distFromSurface <= fgToleranceHalf &&
          facet.Distance(p, fgTolerance) <= fgToleranceHalf)
      {
        // We are on a surface. Return zero.
        aConvex = (fExtremeFacets.find(&facet) != fExtremeFacets.end());
//        Normal(p, aNormalVector);
//        aNormalVector = facet.GetSurfaceNormal();
        aNormalVector = normal;
        return 0.0;
      }
      if (dist >= 0.0 && dist < minDist)
      {
        minDist   = dist;
        minNormal = normal;
        isExtreme = (fExtremeFacets.find(&facet) != fExtremeFacets.end());
        /*
        // FASTER IT IS NOT ...
        bool isExtreme = binary_search (fExtremeFacets2.begin(), fExtremeFacets2.end(), i);
        */
      }
    }
  }
  if (minDist < UUtils::kInfinity)
  {
    aNormalVector = minNormal;
    aConvex = isExtreme;
    return minDist;
  }
  else
  {
    // No intersection found
    aConvex = false;
    Normal(p, aNormalVector);
    return 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::DistanceToOutCandidates(const vector<int>& candidates, const UVector3& aPoint,
                                                const UVector3& direction, double& minDist, UVector3& minNormal, int& minCandidate /*double aPstep,*/ /* UBits & bits*/) const
{
  int candidatesCount = candidates.size();
  double dist            = 0.0;
  double distFromSurface = 0.0;
  UVector3 normal;

  for (int i = 0 ; i < candidatesCount; ++i)
  {
    int candidate = candidates[i];

    VUFacet& facet = *fFacets[candidate];
    if (facet.Intersect(aPoint, direction, true, dist, distFromSurface, normal))
    {
      if (distFromSurface > 0.0 && distFromSurface <= fgToleranceHalf &&
          facet.Distance(aPoint, fgTolerance) <= fgToleranceHalf)
      {
        // We are on a surface
        minDist = 0.0;
        minNormal = normal;
        minCandidate = candidate;
        break;
      }
      if (dist >= 0.0 && dist < minDist)
      {
        minDist = dist;
        minNormal = normal;
        minCandidate = candidate;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::DistanceToOutCore(const UVector3& aPoint, const UVector3& aDirection, UVector3&       aNormalVector, bool& aConvex, double aPstep) const
{
  double minDistance;
  if (fVoxels.GetCountOfVoxels() > 1)
  {
    minDistance = UUtils::kInfinity;

    UVector3 direction = aDirection.Unit();
    double shift = 0;
    vector<int> curVoxel(3);
    if (!fVoxels.Contains(aPoint)) return 0;

    fVoxels.GetVoxel(curVoxel, aPoint);

//    const vector<int> *old = NULL;

//    UBits exclusion (1+0*fVoxels.GetBitsPerSlice());

    int minCandidate = -1;
    do
    {
      const vector<int>& candidates = fVoxels.GetCandidates(curVoxel);
//      if (old == &candidates)
//        old++;

      if (/*old != &candidates &&*/ candidates.size())
      {
        DistanceToOutCandidates(candidates, aPoint, direction, minDistance, aNormalVector, minCandidate);
        if (minDistance <= shift) break;
      }

      shift = fVoxels.DistanceToNext(aPoint, direction, curVoxel);
      if (shift == UUtils::kInfinity) break;

//      old = &candidates;
    }
    while (minDistance > shift);

    if (minCandidate < 0)
    {
      // No intersection found
      minDistance = 0;
      aConvex = false;
      Normal(aPoint, aNormalVector);
    }
    else aConvex = (fExtremeFacets.find(fFacets[minCandidate]) != fExtremeFacets.end());
  }
  else
  {
    minDistance = DistanceToOutNoVoxels(aPoint, aDirection, aNormalVector, aConvex, aPstep);
  }
  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::DistanceToInCandidates(const std::vector<int>& candidates, const UVector3& aPoint, const UVector3& direction /*, double aPstep, UBits & bits*/) const
{
  int candidatesCount = candidates.size();
  double dist            = 0.0;
  double distFromSurface = 0.0;
  UVector3 normal;

  double minDistance = UUtils::kInfinity;
  for (int i = 0 ; i < candidatesCount; ++i)
  {
    int candidate = candidates[i];
//    bits.SetBitNumber(candidate);
    VUFacet& facet = *fFacets[candidate];
    if (facet.Intersect(aPoint, direction, false, dist, distFromSurface, normal))
    {
      //
      // Set minDist to the new distance to current facet if distFromSurface is in
      // positive direction and point is not at surface.  If the point is within
      // 0.5*fgTolerance of the surface, then force distance to be zero and
      // leave member function immediately (for efficiency), as proposed by & credit
      // to Akira Okumura.
      //
      if (distFromSurface > fgToleranceHalf && dist >= 0.0 && dist < minDistance)
      {
        minDistance  = dist;
      }
      if (-fgToleranceHalf <= dist && dist <= fgToleranceHalf)
      {
        return 0.0;
      }
      else
      {
        if (distFromSurface > - fgToleranceHalf && distFromSurface < fgToleranceHalf) minDistance = dist;

      }

    }
  }
  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::DistanceToInCore(const UVector3& aPoint, const UVector3& aDirection, double aPstep) const
{
  double minDistance, distance;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    minDistance = UUtils::kInfinity;
    UVector3 currentPoint = aPoint;
    UVector3 direction = aDirection.Unit();
    double shift = fVoxels.DistanceToFirst(currentPoint, direction);
    if (shift == UUtils::kInfinity) return shift;
    if (shift) currentPoint += direction * shift;
    //    if (!fVoxels.Contains(currentPoint))
    //      return minDistance;

//    UBits exclusion; // (1/*fVoxels.GetBitsPerSlice()*/);
    vector<int> curVoxel(3);

    fVoxels.GetVoxel(curVoxel, currentPoint);
    do
    {
      const vector<int>& candidates = fVoxels.GetCandidates(curVoxel);
      if (candidates.size())
      {
        distance = DistanceToInCandidates(candidates, aPoint, direction);
        if (minDistance > distance) minDistance = distance;
        if (distance < shift) break;
      }
      shift = fVoxels.DistanceToNext(aPoint, direction, curVoxel);
    }
    while (minDistance > shift);

#ifdef DEBUG
    if (fabs(minDistance - distanceToInNoVoxels) > VUSolid::Tolerance())
    {
      VUSolid::EnumInside location = Inside(aPoint);
      minDistance = distanceToInNoVoxels; // you can place a breakpoint here
    }
#endif
//    if (minDistance != UUtils::kInfinity) minDistance += shift;
  }
  else
  {
    minDistance = DistanceToInNoVoxels(aPoint, aDirection, aPstep);
  }

  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
bool UTessellatedSolid::CompareSortedVoxel(const std::pair<int, double>& l,
                                           const std::pair<int, double>& r)
{
  return l.second < r.second;
}

///////////////////////////////////////////////////////////////////////////////
//
// double DistanceToIn(const UVector3& p)
//
// Calculate distance to nearest surface of shape from an outside point p. The
// distance can be an underestimate.
//
double UTessellatedSolid::MinDistanceFacet(const UVector3& p, bool simple, VUFacet*& minFacet) const
{
  double minDist = UUtils::kInfinity;

  int size = fVoxels.GetVoxelBoxesSize();
  vector<pair<int, double> > voxelsSorted(size);

  pair<int, double> info;

  for (int i = 0; i < size; ++i)
  {
    const UVoxelBox& voxelBox = fVoxels.GetVoxelBox(i);

    UVector3 pointShifted = p - voxelBox.pos;
    double safety = fVoxels.MinDistanceToBox(pointShifted, voxelBox.hlen);
    info.first = i;
    info.second = safety;
    voxelsSorted[i] = info;
  }

  std::sort(voxelsSorted.begin(), voxelsSorted.end(), &UTessellatedSolid::CompareSortedVoxel);

  for (int i = 0; i < size; ++i)
  {
    const pair<int, double>& inf = voxelsSorted[i];
//    const UVoxelBox &voxelBox = fVoxels.fVoxelBoxes[inf.first];
    double dist = inf.second;
    if (dist > minDist) break;

    const vector<int>& candidates = fVoxels.GetVoxelBoxCandidates(inf.first);
    int csize = candidates.size();
    for (int j = 0; j < csize; ++j)
    {
      int candidate = candidates[j];
      VUFacet& facet = *fFacets[candidate];
      dist = simple ? facet.Distance(p, minDist) : facet.Distance(p, minDist, false);
      if (dist < minDist)
      {
        minDist  = dist;
        minFacet = &facet;
      }
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::SafetyFromOutside(const UVector3& p, bool aAccurate) const
{
#if UDEBUG
  if (Inside(p) == kInside)
  {
    std::ostringstream message;
    int oldprc = message.precision(16) ;
    message << "Point p is already inside!?" << endl
            << "Position:"  << endl << endl
            << "p.x() = "   << p.x() / mm << " mm" << endl
            << "p.y() = "   << p.y() / mm << " mm" << endl
            << "p.z() = "   << p.z() / mm << " mm" << endl
            << "DistanceToOut(p) == " << DistanceToOut(p);
    message.precision(oldprc) ;
    UUtils::Exception("UTriangularFacet::DistanceToIn(p)", "GeomSolids1002",
                      UWarning, 1, message.str().c_str());
  }
#endif

  double minDist;

  if (!aAccurate) return fVoxels.SafetyToBoundingBox(p);

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    if (!OutsideOfExtent(p, fgTolerance))
    {
      vector<int> startingVoxel(3);
      fVoxels.GetVoxel(startingVoxel, p);
      const vector<int>& candidates = fVoxels.GetCandidates(startingVoxel);
//      int limit = fVoxels.GetCandidatesVoxelArray(p, candidates, NULL);
      if (candidates.size() == 0 && fInsides.GetNbits())
      {
//        int index = fVoxels.GetPointIndex(p);
        int index = fVoxels.GetVoxelsIndex(startingVoxel);
        if (fInsides[index]) return 0.;
      }
    }

    VUFacet* facet;
    minDist = MinDistanceFacet(p, true, facet);
  }
  else
  {
    minDist = UUtils::kInfinity;
    int size = fFacets.size();
    for (int i = 0; i < size; ++i)
    {
      VUFacet& facet = *fFacets[i];
      double dist = facet.Distance(p, minDist);
      if (dist < minDist) minDist  = dist;
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
// double SafetyFrominside(const UVector3& p)
//
// Calculate distance to nearest surface of shape from an inside
// point. The distance can be an underestimate.
//
double UTessellatedSolid::SafetyFromInside(const UVector3& p, bool) const
{
#if UDEBUG
  if (Inside(p) == kOutside)
  {
    std::ostringstream message;
    int oldprc = message.precision(16) ;
    message << "Point p is already outside!?" << endl
            << "Position:"  << endl << endl
            << "p.x() = "   << p.x() / mm << " mm" << endl
            << "p.y() = "   << p.y() / mm << " mm" << endl
            << "p.z() = "   << p.z() / mm << " mm" << endl
            << "DistanceToIn(p) == " << DistanceToIn(p);
    message.precision(oldprc) ;
    UUtils::Exception("UTriangularFacet::DistanceToOut(p)", "GeomSolids1002",
                      UWarning, 1, message.str().c_str());
  }
#endif

  double minDist;

  if (OutsideOfExtent(p, fgTolerance)) return 0.0;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    VUFacet* facet;
    minDist = MinDistanceFacet(p, true, facet);
  }
  else
  {
    minDist = UUtils::kInfinity;
    double dist = 0.0;
    int size = fFacets.size();
    for (int i = 0; i < size; ++i)
    {
      VUFacet& facet = *fFacets[i];
      dist = facet.Distance(p, minDist);
      if (dist < minDist) minDist  = dist;
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
// UGeometryType GetEntityType() const;
//
// Provide identification of the class of an object (required for
// persistency and STEP interface).
//
UGeometryType UTessellatedSolid::GetEntityType() const
{
  return fGeometryType;
}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream& UTessellatedSolid::StreamInfo(std::ostream& os) const
{
  os << endl;
  os << "Geometry Type    = " << fGeometryType  << endl;
  os << "Number of facets = " << fFacets.size() << endl;

  int size = fFacets.size();
  for (int i = 0; i < size; ++i)
  {
    os << "FACET #          = " << i + 1 << endl;
    VUFacet& facet = *fFacets[i];
    facet.StreamInfo(os);
  }
  os << endl;

  return os;
}

///////////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
VUSolid* UTessellatedSolid::Clone() const
{
  return new UTessellatedSolid(*this);
}

///////////////////////////////////////////////////////////////////////////////
//
void UTessellatedSolid::Extent(UVector3& aMin, UVector3& aMax) const
{
  aMin = fMinExtent;
  aMax = fMaxExtent;
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetMinXExtent() const
{
  return fMinExtent.x();
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetMaxXExtent() const
{
  return fMaxExtent.x();
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetMinYExtent() const
{
  return fMinExtent.y();
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetMaxYExtent() const
{
  return fMaxExtent.y();
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetMinZExtent() const
{
  return fMinExtent.z();
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetMaxZExtent() const
{
  return fMaxExtent.z();
}

///////////////////////////////////////////////////////////////////////////////
//
double UTessellatedSolid::GetSurfaceArea()
{
  if (fSurfaceArea != 0.) return fSurfaceArea;

  int size = fFacets.size();
  for (int i = 0; i < size; ++i)
  {
    VUFacet& facet = *fFacets[i];
    fSurfaceArea += facet.GetArea();
  }
  return fSurfaceArea;
}

///////////////////////////////////////////////////////////////////////////////
//
UVector3 UTessellatedSolid::GetPointOnSurface() const
{
  // Select randomly a facet and return a random point on it

  int i = (int) UUtils::Random(0., fFacets.size());
  return fFacets[i]->GetPointOnFace();
}

///////////////////////////////////////////////////////////////////////////////
//
// SetRandomVectorSet
//
// This is a set of predefined random vectors (if that isn't a contradition
// in terms!) used to generate rays from a user-defined point.  The member
// function Inside uses these to determine whether the point is inside or
// outside of the tessellated solid.  All vectors should be unit vectors.
//
void UTessellatedSolid::SetRandomVectors()
{
  fRandir.resize(20);
  fRandir[0] = UVector3(-0.9577428892113370, 0.2732676269591740, 0.0897405271949221);
  fRandir[1]  = UVector3(-0.8331264504940770, -0.5162067214954600, -0.1985722492445700);
  fRandir[2]  = UVector3(-0.1516671651108820, 0.9666292616127460, 0.2064580868390110);
  fRandir[3]  = UVector3(0.6570250350323190, -0.6944539025883300, 0.2933460081893360);
  fRandir[4]  = UVector3(-0.4820456281280320, -0.6331060000098690, -0.6056474264406270);
  fRandir[5]  = UVector3(0.7629032554236800 , 0.1016854697539910, -0.6384658864065180);
  fRandir[6]  = UVector3(0.7689540409061150, 0.5034929891988220, 0.3939600142169160);
  fRandir[7]  = UVector3(0.5765188359255740, 0.5997271636278330, -0.5549354566343150);
  fRandir[8]  = UVector3(0.6660632777862070, -0.6362809868288380, 0.3892379937580790);
  fRandir[9]  = UVector3(0.3824415020414780, 0.6541792713761380, -0.6525243125110690);
  fRandir[10] = UVector3(-0.5107726564526760, 0.6020905056811610, 0.6136760679616570);
  fRandir[11] = UVector3(0.7459135439578050, 0.6618796061649330, 0.0743530220183488);
  fRandir[12] = UVector3(0.1536405855311580, 0.8117477913978260, -0.5634359711967240);
  fRandir[13] = UVector3(0.0744395301705579, -0.8707110101772920, -0.4861286795736560);
  fRandir[14] = UVector3(-0.1665874645185400, 0.6018553940549240, -0.7810369397872780);
  fRandir[15] = UVector3(0.7766902003633100, 0.6014617505959970, -0.1870724331097450);
  fRandir[16] = UVector3(-0.8710128685847430, -0.1434320216603030, -0.4698551243971010);
  fRandir[17] = UVector3(0.8901082092766820, -0.4388411398893870, 0.1229871120030100);
  fRandir[18] = UVector3(-0.6430417431544370, -0.3295938228697690, 0.6912779675984150);
  fRandir[19] = UVector3(0.6331124368380410, 0.6306211461665000, 0.4488714875425340);

  fMaxTries = 20;
}

double const UTessellatedSolid::dirTolerance = 1.0E-14;

///////////////////////////////////////////////////////////////////////////////
//
int UTessellatedSolid::AllocatedMemoryWithoutVoxels()
{
  int base = sizeof(*this);
  base += fVertexList.capacity() * sizeof(UVector3);
  base += fRandir.capacity() * sizeof(UVector3);

  int limit = fFacets.size();
  for (int i = 0; i < limit; ++i)
  {
    VUFacet& facet = *fFacets[i];
    base += facet.AllocatedMemory() + sizeof(VUFacet*);
  }

  std::set<VUFacet*>::const_iterator beg, end, it;
  beg = fExtremeFacets.begin();
  end = fExtremeFacets.end();
  for (it = beg; it != end; it++)
  {
    VUFacet& facet = *(*it);
    base += facet.AllocatedMemory();
  }
  return base;
}

///////////////////////////////////////////////////////////////////////////////
//
int UTessellatedSolid::AllocatedMemory()
{
  int size = AllocatedMemoryWithoutVoxels();
  int sizeInsides = fInsides.GetNbytes();
  int sizeVoxels = fVoxels.AllocatedMemory();
  size += sizeInsides + sizeVoxels;
  return size;
}
