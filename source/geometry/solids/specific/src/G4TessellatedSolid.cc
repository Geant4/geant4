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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TessellatedSolid.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
// 23 October 2016,   E Tcherniaev, reimplemented CalculateExtent() to make
//                    use of G4BoundingEnvelope.
//
// 12 October 2012,   M Gayer, CERN, complete rewrite reducing memory
//                    requirements more than 50% and speedup by a factor of
//                    tens or more depending on the number of facets, thanks
//                    to voxelization of surface and improvements. 
//                    Speedup factor of thousands for solids with number of
//                    facets in hundreds of thousands.
//
// 22 August 2011,    I Hrivnacova, Orsay, fix in DistanceToOut(p) and
//                    DistanceToIn(p) to exactly compute distance from facet
//                    avoiding use of 'outgoing' flag shortcut variant.
//
// 04 August 2011,    T Nikitina, CERN, added SetReferences() to
//                    CreatePolyhedron() for Visualization of Boolean Operations  
//
// 12 April 2010,     P R Truscott, QinetiQ, bug fixes to treat optical
//                    photon transport, in particular internal reflection
//                    at surface.
//
// 14 November 2007,  P R Truscott, QinetiQ & Stan Seibert, U Texas
//                    Bug fixes to CalculateExtent
//
// 17 September 2007, P R Truscott, QinetiQ Ltd & Richard Holmberg
//                    Updated extensively prior to this date to deal with
//                    concaved tessellated surfaces, based on the algorithm
//                    of Richard Holmberg.  This had been slightly modified
//                    to determine with inside the geometry by projecting
//                    random rays from the point provided.  Now random rays
//                    are predefined rather than making use of random
//                    number generator at run-time.
//
// 22 November 2005, F Lei
//  - Changed ::DescribeYourselfTo(), line 464
//  - added GetPolyHedron()
// 
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include <iostream>
#include <stack>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <list>

#include "G4TessellatedSolid.hh"

#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeometryTolerance.hh"
#include "G4VFacet.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4PolyhedronArbitrary.hh"
#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// Standard contructor has blank name and defines no fFacets.
//
G4TessellatedSolid::G4TessellatedSolid () : G4VSolid("dummy")
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor. Simple define name and geometry type - no fFacets
// to detine.
//
G4TessellatedSolid::G4TessellatedSolid (const G4String &name)
  : G4VSolid(name)
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4TessellatedSolid::G4TessellatedSolid( __void__& a) : G4VSolid(a)
{
  Initialize();
  fMinExtent.set(0,0,0);
  fMaxExtent.set(0,0,0);
}


///////////////////////////////////////////////////////////////////////////////
G4TessellatedSolid::~G4TessellatedSolid ()
{
  DeleteObjects ();
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor.
//
G4TessellatedSolid::G4TessellatedSolid (const G4TessellatedSolid &ts)
  : G4VSolid(ts), fpPolyhedron(0)
{
  Initialize();

  CopyObjects(ts);
}

///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator.
//
G4TessellatedSolid&
G4TessellatedSolid::operator= (const G4TessellatedSolid &ts)
{
  if (&ts == this) return *this;

  // Copy base class data
  G4VSolid::operator=(ts);

  DeleteObjects ();

  Initialize();

  CopyObjects (ts);

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::Initialize()
{
  kCarToleranceHalf = 0.5*kCarTolerance;

  fRebuildPolyhedron = false; fpPolyhedron = 0;
  fCubicVolume = 0.; fSurfaceArea = 0.;

  fGeometryType = "G4TessellatedSolid";
  fSolidClosed  = false;

  fMinExtent.set(kInfinity,kInfinity,kInfinity);
  fMaxExtent.set(-kInfinity,-kInfinity,-kInfinity);

  SetRandomVectors();
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DeleteObjects ()
{
  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)  { delete fFacets[i]; }
  fFacets.clear();
  delete fpPolyhedron; fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::CopyObjects (const G4TessellatedSolid &ts)
{
  G4ThreeVector reductionRatio;
  G4int fmaxVoxels = fVoxels.GetMaxVoxels(reductionRatio);
  if (fmaxVoxels < 0)
    fVoxels.SetMaxVoxels(reductionRatio);
  else
    fVoxels.SetMaxVoxels(fmaxVoxels);

  G4int n = ts.GetNumberOfFacets();
  for (G4int i = 0; i < n; ++i)
  {
    G4VFacet *facetClone = (ts.GetFacet(i))->GetClone();
    AddFacet(facetClone);
  }
  if (ts.GetSolidClosed()) SetSolidClosed(true);
}

///////////////////////////////////////////////////////////////////////////////
//
// Add a facet to the facet list.
// Note that you can add, but you cannot delete.
//
G4bool G4TessellatedSolid::AddFacet (G4VFacet *aFacet)
{
  // Add the facet to the vector.
  //
  if (fSolidClosed)
  {
    G4Exception("G4TessellatedSolid::AddFacet()", "GeomSolids1002",
                JustWarning, "Attempt to add facets when solid is closed.");
    return false;
  }
  else if (aFacet->IsDefined())
  {
    set<G4VertexInfo,G4VertexComparator>::iterator begin
      = fFacetList.begin(), end = fFacetList.end(), pos, it;
    G4ThreeVector p = aFacet->GetCircumcentre();
    G4VertexInfo value;
    value.id = fFacetList.size();
    value.mag2 = p.x() + p.y() + p.z();

    G4bool found = false;
    if (!OutsideOfExtent(p, kCarTolerance))
    {
      G4double kCarTolerance3 = 3 * kCarTolerance;
      pos = fFacetList.lower_bound(value);

      it = pos;
      while (!found && it != end)    // Loop checking, 13.08.2015, G.Cosmo
      {
        G4int id = (*it).id;
        G4VFacet *facet = fFacets[id];
        G4ThreeVector q = facet->GetCircumcentre();
        if ((found = (facet == aFacet))) break;
        G4double dif = q.x() + q.y() + q.z() - value.mag2;
        if (dif > kCarTolerance3) break;
        it++;
      }

      if (fFacets.size() > 1)
      {
        it = pos;
        while (!found && it != begin)    // Loop checking, 13.08.2015, G.Cosmo
        {
          --it;
          G4int id = (*it).id;
          G4VFacet *facet = fFacets[id];  
          G4ThreeVector q = facet->GetCircumcentre();
          found = (facet == aFacet);
          if (found) break;
          G4double dif = value.mag2 - (q.x() + q.y() + q.z());
          if (dif > kCarTolerance3) break;
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
    G4Exception("G4TessellatedSolid::AddFacet()", "GeomSolids1002",
                JustWarning, "Attempt to add facet not properly defined.");    
    aFacet->StreamInfo(G4cout);
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::SetAllUsingStack(const std::vector<G4int> &voxel,
                                           const std::vector<G4int> &max,
                                           G4bool status, G4SurfBits &checked)
{
  vector<G4int> xyz = voxel;
  stack<vector<G4int> > pos;
  pos.push(xyz);
  G4int filled = 0;
  G4int cc = 0, nz = 0;

  vector<G4int> candidates;

  while (!pos.empty())    // Loop checking, 13.08.2015, G.Cosmo
  {
    xyz = pos.top();
    pos.pop();
    G4int index = fVoxels.GetVoxelsIndex(xyz);
    if (!checked[index])
    {
      checked.SetBitNumber(index, true);
      cc++;

      if (fVoxels.IsEmpty(index))
      {
        filled++;

        fInsides.SetBitNumber(index, status);

        for (G4int i = 0; i <= 2; ++i)
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
void G4TessellatedSolid::PrecalculateInsides()
{
  vector<G4int> voxel(3), maxVoxels(3);
  for (G4int i = 0; i <= 2; ++i) maxVoxels[i] = fVoxels.GetBoundary(i).size();
  G4int size = maxVoxels[0] * maxVoxels[1] * maxVoxels[2];

  G4SurfBits checked(size-1);
  fInsides.Clear();
  fInsides.ResetBitNumber(size-1);

  G4ThreeVector point;
  for (voxel[2] = 0; voxel[2] < maxVoxels[2] - 1; ++voxel[2])
  {
    for (voxel[1] = 0; voxel[1] < maxVoxels[1] - 1; ++voxel[1])
    {
      for (voxel[0] = 0; voxel[0] < maxVoxels[0] - 1; ++voxel[0])
      {
        G4int index = fVoxels.GetVoxelsIndex(voxel);
        if (!checked[index] && fVoxels.IsEmpty(index))
        {
          for (G4int i = 0; i <= 2; ++i)
          {
            point[i] = fVoxels.GetBoundary(i)[voxel[i]];
          }
          G4bool inside = (G4bool) (InsideNoVoxels(point) == kInside);
          SetAllUsingStack(voxel, maxVoxels, inside, checked);
        }
        else checked.SetBitNumber(index);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::Voxelize ()
{
#ifdef G4SPECSDEBUG
  G4cout << "Voxelizing..." << G4endl;
#endif
  fVoxels.Voxelize(fFacets);

  if (fVoxels.Empty().GetNbits())
  {
#ifdef G4SPECSDEBUG
    G4cout << "Precalculating Insides..." << G4endl;
#endif
    PrecalculateInsides();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Compute extremeFacets, i.e. find those facets that have surface
// planes that bound the volume.
// Note that this is going to reject concaved surfaces as being extreme. Also
// note that if the vertex is on the facet, displacement is zero, so IsInside
// returns true. So will this work??  Need non-equality
// "G4bool inside = displacement < 0.0;"
// or
// "G4bool inside = displacement <= -0.5*kCarTolerance" 
// (Notes from PT 13/08/2007).
//
void G4TessellatedSolid::SetExtremeFacets()
{
  G4int size = fFacets.size();
  for (G4int j = 0; j < size; ++j)
  {
    G4VFacet &facet = *fFacets[j];

    G4bool isExtreme = true;
    G4int vsize = fVertexList.size();
    for (G4int i=0; i < vsize; ++i)
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
void G4TessellatedSolid::CreateVertexList()
{
  // The algorithm:
  // we will have additional vertexListSorted, where all the items will be
  // sorted by magnitude of vertice vector.
  // New candidate for fVertexList - we will determine the position fo first
  // item which would be within its magnitude - 0.5*kCarTolerance.
  // We will go trough until we will reach > +0.5 kCarTolerance.
  // Comparison (q-p).mag() < 0.5*kCarTolerance will be made.
  // They can be just stored in std::vector, with custom insertion based
  // on binary search.

  set<G4VertexInfo,G4VertexComparator> vertexListSorted;
  set<G4VertexInfo,G4VertexComparator>::iterator begin
     = vertexListSorted.begin(), end = vertexListSorted.end(), pos, it;
  G4ThreeVector p;
  G4VertexInfo value;

  fVertexList.clear();
  G4int size = fFacets.size();

  G4double kCarTolerance24 = kCarTolerance * kCarTolerance / 4.0;
  G4double kCarTolerance3 = 3 * kCarTolerance;
  vector<G4int> newIndex(100);
  
  for (G4int k = 0; k < size; ++k)
  {
    G4VFacet &facet = *fFacets[k];
    G4int max = facet.GetNumberOfVertices();

    for (G4int i = 0; i < max; ++i)
    {
      p = facet.GetVertex(i);
      value.id = fVertexList.size();
      value.mag2 = p.x() + p.y() + p.z();

      G4bool found = false;
      G4int id = 0;
      if (!OutsideOfExtent(p, kCarTolerance))
      {
        pos = vertexListSorted.lower_bound(value);
        it = pos;
        while (it != end)    // Loop checking, 13.08.2015, G.Cosmo
        {
          id = (*it).id;
          G4ThreeVector q = fVertexList[id];
          G4double dif = (q-p).mag2();
          found = (dif < kCarTolerance24);
          if (found) break;
          dif = q.x() + q.y() + q.z() - value.mag2;
          if (dif > kCarTolerance3) break;
          it++;
        }

        if (!found && (fVertexList.size() > 1))
        {
          it = pos;
          while (it != begin)    // Loop checking, 13.08.2015, G.Cosmo
          {
            --it;
            id = (*it).id;
            G4ThreeVector q = fVertexList[id];
            G4double dif = (q-p).mag2();
            found = (dif < kCarTolerance24);
            if (found) break;
        dif = value.mag2 - (q.x() + q.y() + q.z());
            if (dif > kCarTolerance3) break;
          }
        }
      }

      if (!found)
      {
#ifdef G4SPECSDEBUG
        G4cout << p.x() << ":" << p.y() << ":" << p.z() << G4endl;
        G4cout << "Adding new vertex #" << i << " of facet " << k
               << " id " << value.id << G4endl;
        G4cout << "===" << G4endl;
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
          if (p.x() > fMaxExtent.x()) fMaxExtent.setX(p.x());
          else if (p.x() < fMinExtent.x()) fMinExtent.setX(p.x());
          if (p.y() > fMaxExtent.y()) fMaxExtent.setY(p.y());
          else if (p.y() < fMinExtent.y()) fMinExtent.setY(p.y());
          if (p.z() > fMaxExtent.z()) fMaxExtent.setZ(p.z());
          else if (p.z() < fMinExtent.z()) fMinExtent.setZ(p.z());
        }
      }
      else
      {
#ifdef G4SPECSDEBUG
        G4cout << p.x() << ":" << p.y() << ":" << p.z() << G4endl;
        G4cout << "Vertex #" << i << " of facet " << k
               << " found, redirecting to " << id << G4endl;
        G4cout << "===" << G4endl;
#endif
        newIndex[i] = id;
      }
    }
    // only now it is possible to change vertices pointer
    //
    facet.SetVertices(&fVertexList);
    for (G4int i = 0; i < max; i++)
        facet.SetVertexIndex(i,newIndex[i]);
  }
  vector<G4ThreeVector>(fVertexList).swap(fVertexList);
  
#ifdef G4SPECSDEBUG
  G4double previousValue = 0;
  for (set<G4VertexInfo,G4VertexComparator>::iterator res=
       vertexListSorted.begin(); res!=vertexListSorted.end(); ++res)
  {
    G4int id = (*res).id;
    G4ThreeVector vec = fVertexList[id];
    G4double mvalue = vec.x() + vec.y() + vec.z();
    if (previousValue && (previousValue - 1e-9 > mvalue))
      G4cout << "Error in CreateVertexList: previousValue " << previousValue 
             <<  " is smaller than mvalue " << mvalue << G4endl;
    previousValue = mvalue;
  }
#endif
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DisplayAllocatedMemory()
{
  G4int without = AllocatedMemoryWithoutVoxels();
  G4int with = AllocatedMemory();
  G4double ratio = (G4double) with / without;
  G4cout << "G4TessellatedSolid - Allocated memory without voxel overhead "
         << without << "; with " << with << "; ratio: " << ratio << G4endl; 
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::SetSolidClosed (const G4bool t)
{
  if (t)
  {
#ifdef G4SPECSDEBUG    
    G4cout << "Creating vertex list..." << G4endl;
#endif
    CreateVertexList();

#ifdef G4SPECSDEBUG    
    G4cout << "Setting extreme facets..." << G4endl;
#endif
    SetExtremeFacets();
    
#ifdef G4SPECSDEBUG    
    G4cout << "Voxelizing..." << G4endl;
#endif
    Voxelize();

#ifdef G4SPECSDEBUG
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
G4bool G4TessellatedSolid::GetSolidClosed () const
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
G4TessellatedSolid &
G4TessellatedSolid::operator+=(const G4TessellatedSolid &right)
{
  G4int size = right.GetNumberOfFacets();
  for (G4int i = 0; i < size; ++i)
    AddFacet(right.GetFacet(i)->GetClone());

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetNumberOfFacets
//
G4int G4TessellatedSolid::GetNumberOfFacets () const
{
  return fFacets.size();
}

///////////////////////////////////////////////////////////////////////////////
//
EInside G4TessellatedSolid::InsideVoxels(const G4ThreeVector &p) const
{
  //
  // First the simple test - check if we're outside of the X-Y-Z extremes
  // of the tessellated solid.
  //
  if (OutsideOfExtent(p, kCarTolerance))
    return kOutside;

  vector<G4int> startingVoxel(3);
  fVoxels.GetVoxel(startingVoxel, p);

  const G4double dirTolerance = 1.0E-14;

  const vector<G4int> &startingCandidates =
    fVoxels.GetCandidates(startingVoxel);
  G4int limit = startingCandidates.size();
  if (limit == 0 && fInsides.GetNbits())
  {
    G4int index = fVoxels.GetPointIndex(p);
    EInside location = fInsides[index] ? kInside : kOutside;
    return location;
  }

  G4double minDist = kInfinity;

  for(G4int i = 0; i < limit; ++i)
  {
    G4int candidate = startingCandidates[i];
    G4VFacet &facet = *fFacets[candidate];
    G4double dist = facet.Distance(p,minDist);
    if (dist < minDist) minDist = dist;
    if (dist <= kCarToleranceHalf)
      return kSurface;
  }

  // The following is something of an adaptation of the method implemented by
  // Rickard Holmberg augmented with information from Schneider & Eberly,
  // "Geometric Tools for Computer Graphics," pp700-701, 2003. In essence,
  // we're trying to determine whether we're inside the volume by projecting
  // a few rays and determining if the first surface crossed is has a normal
  // vector between 0 to pi/2 (out-going) or pi/2 to pi (in-going).
  // We should also avoid rays which are nearly within the plane of the
  // tessellated surface, and therefore produce rays randomly.
  // For the moment, this is a bit over-engineered (belt-braces-and-ducttape).
  //
  G4double distOut          = kInfinity;
  G4double distIn           = kInfinity;
  G4double distO            = 0.0;
  G4double distI            = 0.0;
  G4double distFromSurfaceO = 0.0;
  G4double distFromSurfaceI = 0.0;
  G4ThreeVector normalO, normalI;
  G4bool crossingO          = false;
  G4bool crossingI          = false;
  EInside location          = kOutside;
  G4int sm                  = 0;

  G4bool nearParallel = false;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    // We loop until we find direction where the vector is not nearly parallel
    // to the surface of any facet since this causes ambiguities.  The usual
    // case is that the angles should be sufficiently different, but there
    // are 20 random directions to select from - hopefully sufficient.
    //
    distOut = distIn = kInfinity;
    const G4ThreeVector &v = fRandir[sm];
    sm++;
    //
    // This code could be voxelized by the same algorithm, which is used for
    // DistanceToOut(). We will traverse through fVoxels. we will call
    // intersect only for those, which would be candidates and was not
    // checked before.
    //
    G4ThreeVector currentPoint = p;
    G4ThreeVector direction = v.unit();
    // G4SurfBits exclusion(fVoxels.GetBitsPerSlice());
    vector<G4int> curVoxel(3);
    curVoxel = startingVoxel;
    G4double shiftBonus = kCarTolerance;

    G4bool crossed = false;
    G4bool started = true;

    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      const vector<G4int> &candidates =
        started ? startingCandidates : fVoxels.GetCandidates(curVoxel);
      started = false;
      if (G4int candidatesCount = candidates.size())
      {  
        for (G4int i = 0 ; i < candidatesCount; ++i)
        {
          G4int candidate = candidates[i];
          // bits.SetBitNumber(candidate);
          G4VFacet &facet = *fFacets[candidate];

          crossingO = facet.Intersect(p,v,true,distO,distFromSurfaceO,normalO);
          crossingI = facet.Intersect(p,v,false,distI,distFromSurfaceI,normalI);

          if (crossingO || crossingI)
          {
            crossed = true;

            nearParallel = (crossingO
                     && std::fabs(normalO.dot(v))<dirTolerance)
                     || (crossingI && std::fabs(normalI.dot(v))<dirTolerance);
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
          G4int index = fVoxels.GetVoxelsIndex(curVoxel);
          G4bool inside = fInsides[index];
          location = inside ? kInside : kOutside;
          return location;
        }
      }

      G4double shift=fVoxels.DistanceToNext(currentPoint, direction, curVoxel);
      if (shift == kInfinity) break;

      currentPoint += direction * (shift + shiftBonus);
    }
    while (fVoxels.UpdateCurrentVoxel(currentPoint, direction, curVoxel));

  }
  while (nearParallel && sm!=fMaxTries);
  //
  // Here we loop through the facets to find out if there is an intersection
  // between the ray and that facet.  The test if performed separately whether
  // the ray is entering the facet or exiting.
  //
#ifdef G4VERBOSE
  if (sm == fMaxTries)
  {
    //
    // We've run out of random vector directions. If nTries is set sufficiently
    // low (nTries <= 0.5*maxTries) then this would indicate that there is
    // something wrong with geometry.
    //
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Cannot determine whether point is inside or outside volume!"
      << G4endl
      << "Solid name       = " << GetName()  << G4endl
      << "Geometry Type    = " << fGeometryType  << G4endl
      << "Number of facets = " << fFacets.size() << G4endl
      << "Position:"  << G4endl << G4endl
      << "p.x() = "   << p.x()/mm << " mm" << G4endl
      << "p.y() = "   << p.y()/mm << " mm" << G4endl
      << "p.z() = "   << p.z()/mm << " mm";
    message.precision(oldprc);
    G4Exception("G4TessellatedSolid::Inside()",
                "GeomSolids1002", JustWarning, message);
  }
#endif

  // In the next if-then-elseif G4String the logic is as follows:
  // (1) You don't hit anything so cannot be inside volume, provided volume
  //     constructed correctly!
  // (2) Distance to inside (ie. nearest facet such that you enter facet) is
  //     shorter than distance to outside (nearest facet such that you exit
  //     facet) - on condition of safety distance - therefore we're outside.
  // (3) Distance to outside is shorter than distance to inside therefore
  //     we're inside.
  //
  if (distIn == kInfinity && distOut == kInfinity)
    location = kOutside;
  else if (distIn <= distOut - kCarToleranceHalf)
    location = kOutside;
  else if (distOut <= distIn - kCarToleranceHalf)
    location = kInside;

  return location;
}
 
///////////////////////////////////////////////////////////////////////////////
//
EInside G4TessellatedSolid::InsideNoVoxels (const G4ThreeVector &p) const
{
  //
  // First the simple test - check if we're outside of the X-Y-Z extremes
  // of the tessellated solid.
  //
  if (OutsideOfExtent(p, kCarTolerance))
    return kOutside;

  const G4double dirTolerance = 1.0E-14;

  G4double minDist = kInfinity;
  //
  // Check if we are close to a surface
  //
  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    G4double dist = facet.Distance(p,minDist);
    if (dist < minDist) minDist = dist;
    if (dist <= kCarToleranceHalf)
    {
      return kSurface;
    }
  }
  //
  // The following is something of an adaptation of the method implemented by
  // Rickard Holmberg augmented with information from Schneider & Eberly,
  // "Geometric Tools for Computer Graphics," pp700-701, 2003. In essence, we're
  // trying to determine whether we're inside the volume by projecting a few
  // rays and determining if the first surface crossed is has a normal vector
  // between 0 to pi/2 (out-going) or pi/2 to pi (in-going). We should also
  // avoid rays which are nearly within the plane of the tessellated surface,
  // and therefore produce rays randomly. For the moment, this is a bit
  // over-engineered (belt-braces-and-ducttape).
  //
#if G4SPECSDEBUG
  G4int nTry                = 7;
#else
  G4int nTry                = 3;
#endif
  G4double distOut          = kInfinity;
  G4double distIn           = kInfinity;
  G4double distO            = 0.0;
  G4double distI            = 0.0;
  G4double distFromSurfaceO = 0.0;
  G4double distFromSurfaceI = 0.0;
  G4ThreeVector normalO(0.0,0.0,0.0);
  G4ThreeVector normalI(0.0,0.0,0.0);
  G4bool crossingO          = false;
  G4bool crossingI          = false;
  EInside location          = kOutside;
  EInside locationprime     = kOutside;
  G4int sm = 0;

  for (G4int i=0; i<nTry; ++i)
  {
    G4bool nearParallel = false;
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      //
      // We loop until we find direction where the vector is not nearly parallel
      // to the surface of any facet since this causes ambiguities.  The usual
      // case is that the angles should be sufficiently different, but there
      // are 20 random directions to select from - hopefully sufficient.
      //
      distOut =  distIn = kInfinity;
      G4ThreeVector v = fRandir[sm];
      sm++;
      vector<G4VFacet*>::const_iterator f = fFacets.begin();

      do    // Loop checking, 13.08.2015, G.Cosmo
      {
        //
        // Here we loop through the facets to find out if there is an
        // intersection between the ray and that facet. The test if performed
        // separately whether the ray is entering the facet or exiting.
        //
        crossingO = ((*f)->Intersect(p,v,true,distO,distFromSurfaceO,normalO));
        crossingI = ((*f)->Intersect(p,v,false,distI,distFromSurfaceI,normalI));
        if (crossingO || crossingI)
        {
          nearParallel = (crossingO && std::fabs(normalO.dot(v))<dirTolerance)
                      || (crossingI && std::fabs(normalI.dot(v))<dirTolerance);
          if (!nearParallel)
          {
            if (crossingO && distO > 0.0 && distO < distOut) distOut = distO;
            if (crossingI && distI > 0.0 && distI < distIn)  distIn  = distI;
          }
        }
      } while (!nearParallel && ++f!=fFacets.end());
    } while (nearParallel && sm!=fMaxTries);

#ifdef G4VERBOSE
    if (sm == fMaxTries)
    {
      //
      // We've run out of random vector directions. If nTries is set
      // sufficiently low (nTries <= 0.5*maxTries) then this would indicate
      // that there is something wrong with geometry.
      //
      std::ostringstream message;
      G4int oldprc = message.precision(16);
      message << "Cannot determine whether point is inside or outside volume!"
        << G4endl
        << "Solid name       = " << GetName()  << G4endl
        << "Geometry Type    = " << fGeometryType  << G4endl
        << "Number of facets = " << fFacets.size() << G4endl
        << "Position:"  << G4endl << G4endl
        << "p.x() = "   << p.x()/mm << " mm" << G4endl
        << "p.y() = "   << p.y()/mm << " mm" << G4endl
        << "p.z() = "   << p.z()/mm << " mm";
      message.precision(oldprc);
      G4Exception("G4TessellatedSolid::Inside()",
        "GeomSolids1002", JustWarning, message);
    }
#endif
    //
    // In the next if-then-elseif G4String the logic is as follows:
    // (1) You don't hit anything so cannot be inside volume, provided volume
    //     constructed correctly!
    // (2) Distance to inside (ie. nearest facet such that you enter facet) is
    //     shorter than distance to outside (nearest facet such that you exit
    //     facet) - on condition of safety distance - therefore we're outside.
    // (3) Distance to outside is shorter than distance to inside therefore
    // we're inside.
    //
    if (distIn == kInfinity && distOut == kInfinity)
      locationprime = kOutside;
    else if (distIn <= distOut - kCarToleranceHalf)
      locationprime = kOutside;
    else if (distOut <= distIn - kCarToleranceHalf)
      locationprime = kInside;

    if (i == 0) location = locationprime;
  }

  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.
//
G4bool G4TessellatedSolid::Normal (const G4ThreeVector &p,
                                         G4ThreeVector &aNormal) const
{
  G4double minDist;
  G4VFacet *facet = 0;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    vector<G4int> curVoxel(3);
    fVoxels.GetVoxel(curVoxel, p);
    const vector<G4int> &candidates = fVoxels.GetCandidates(curVoxel);
    // fVoxels.GetCandidatesVoxelArray(p, candidates, 0);

    if (G4int limit = candidates.size())
    {
      minDist = kInfinity;
      for(G4int i = 0 ; i < limit ; ++i)
      {      
        G4int candidate = candidates[i];
        G4VFacet &fct = *fFacets[candidate];
        G4double dist = fct.Distance(p,minDist);
        if (dist < minDist) minDist = dist;
        if (dist <= kCarToleranceHalf)
        {
          aNormal = fct.GetSurfaceNormal();
          return true;
        }
      }
    }
    minDist = MinDistanceFacet(p, true, facet);
  }
  else
  {
    minDist = kInfinity;
    G4int size = fFacets.size();
    for (G4int i = 0; i < size; ++i)
    {
      G4VFacet &f = *fFacets[i];
      G4double dist = f.Distance(p, minDist);
      if (dist < minDist)
      {
        minDist  = dist;
        facet = &f;
      }
    }
  }

  if (minDist != kInfinity)
  {
    if (facet)  { aNormal = facet->GetSurfaceNormal(); }
    return minDist <= kCarToleranceHalf;
  }
  else
  {
#ifdef G4VERBOSE
    std::ostringstream message;
    message << "Point p is not on surface !?" << G4endl
      << "          No facets found for point: " << p << " !" << G4endl
      << "          Returning approximated value for normal.";

    G4Exception("G4TessellatedSolid::SurfaceNormal(p)",
                "GeomSolids1002", JustWarning, message );
#endif
    aNormal = (p.z() > 0 ? G4ThreeVector(0,0,1) : G4ThreeVector(0,0,-1));
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v)
//
// Return the distance along the normalised vector v to the shape,
// from the point at offset p. If there is no intersection, return
// kInfinity. The first intersection resulting from 'leaving' a
// surface/volume is discarded. Hence, this is tolerant of points on
// surface of shape.
//
G4double
G4TessellatedSolid::DistanceToInNoVoxels (const G4ThreeVector &p,
                                          const G4ThreeVector &v,
                                                G4double /*aPstep*/) const
{
  G4double minDist         = kInfinity;
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal;

#if G4SPECSDEBUG
  if (Inside(p) == kInside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16) ;
    message << "Point p is already inside!?" << G4endl
      << "Position:"  << G4endl << G4endl
      << "   p.x() = "   << p.x()/mm << " mm" << G4endl
      << "   p.y() = "   << p.y()/mm << " mm" << G4endl
      << "   p.z() = "   << p.z()/mm << " mm" << G4endl
      << "DistanceToOut(p) == " << DistanceToOut(p);
    message.precision(oldprc) ;
    G4Exception("G4TriangularFacet::DistanceToIn(p,v)",
                "GeomSolids1002", JustWarning, message);
  }
#endif

  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    if (facet.Intersect(p,v,false,dist,distFromSurface,normal))
    {
      //
      // set minDist to the new distance to current facet if distFromSurface
      // is in positive direction and point is not at surface. If the point is
      // within 0.5*kCarTolerance of the surface, then force distance to be
      // zero and leave member function immediately (for efficiency), as
      // proposed by & credit to Akira Okumura.
      //
      if (distFromSurface > kCarToleranceHalf && dist >= 0.0 && dist < minDist)
      {
        minDist  = dist;
      }
      else
      {
        if (-kCarToleranceHalf <= dist && dist <= kCarToleranceHalf)
        {
          return 0.0;
        }
        else
        {
          if  (distFromSurface > -kCarToleranceHalf
            && distFromSurface <  kCarToleranceHalf)
          {
            minDist = dist;
          }
        }
      }
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double
G4TessellatedSolid::DistanceToOutNoVoxels (const G4ThreeVector &p,
                                           const G4ThreeVector &v,
                                                 G4ThreeVector &aNormalVector,
                                                 G4bool &aConvex,
                                                 G4double /*aPstep*/) const
{
  G4double minDist         = kInfinity;
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal, minNormal;

#if G4SPECSDEBUG
  if ( Inside(p) == kOutside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16) ;
    message << "Point p is already outside!?" << G4endl
      << "Position:"  << G4endl << G4endl
      << "   p.x() = "   << p.x()/mm << " mm" << G4endl
      << "   p.y() = "   << p.y()/mm << " mm" << G4endl
      << "   p.z() = "   << p.z()/mm << " mm" << G4endl
      << "DistanceToIn(p) == " << DistanceToIn(p);
    message.precision(oldprc) ;
    G4Exception("G4TriangularFacet::DistanceToOut(p)",
                "GeomSolids1002", JustWarning, message);
  }
#endif

  G4bool isExtreme = false;
  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    if (facet.Intersect(p,v,true,dist,distFromSurface,normal))
    {
      if (distFromSurface > 0.0 && distFromSurface <= kCarToleranceHalf &&
        facet.Distance(p,kCarTolerance) <= kCarToleranceHalf)
      {
        // We are on a surface. Return zero.
        aConvex = (fExtremeFacets.find(&facet) != fExtremeFacets.end());
        // Normal(p, aNormalVector);
        // aNormalVector = facet.GetSurfaceNormal();
        aNormalVector = normal;
        return 0.0;
      }
      if (dist >= 0.0 && dist < minDist)
      {
        minDist   = dist;
        minNormal = normal;
        isExtreme = (fExtremeFacets.find(&facet) != fExtremeFacets.end());
      }
    }
  }
  if (minDist < kInfinity)
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
void G4TessellatedSolid::
DistanceToOutCandidates(const std::vector<G4int> &candidates,
                        const G4ThreeVector &aPoint,
                        const G4ThreeVector &direction,
                              G4double &minDist, G4ThreeVector &minNormal,
                              G4int &minCandidate ) const
{
  G4int candidatesCount = candidates.size();
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal;

  for (G4int i = 0 ; i < candidatesCount; ++i)
  {
    G4int candidate = candidates[i];
    G4VFacet &facet = *fFacets[candidate];
    if (facet.Intersect(aPoint,direction,true,dist,distFromSurface,normal))
    {
      if (distFromSurface > 0.0 && distFromSurface <= kCarToleranceHalf
       && facet.Distance(aPoint,kCarTolerance) <= kCarToleranceHalf)
      {
        // We are on a surface
        //
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
G4double
G4TessellatedSolid::DistanceToOutCore(const G4ThreeVector &aPoint,
                                      const G4ThreeVector &aDirection,
                                            G4ThreeVector &aNormalVector,
                                            G4bool &aConvex,
                                            G4double aPstep) const
{
  G4double minDistance;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    minDistance = kInfinity;

    G4ThreeVector currentPoint = aPoint;
    G4ThreeVector direction = aDirection.unit();
    G4double totalShift = 0;
    vector<G4int> curVoxel(3);
    if (!fVoxels.Contains(aPoint)) return 0;

    fVoxels.GetVoxel(curVoxel, currentPoint);

    G4double shiftBonus = kCarTolerance;

    const vector<G4int> *old = 0;

    G4int minCandidate = -1;
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      const vector<G4int> &candidates = fVoxels.GetCandidates(curVoxel);
      if (old == &candidates)
        old++;
      if (old != &candidates && candidates.size())
      {
        DistanceToOutCandidates(candidates, aPoint, direction, minDistance,
                                aNormalVector, minCandidate); 
        if (minDistance <= totalShift) break; 
      }

      G4double shift=fVoxels.DistanceToNext(currentPoint, direction, curVoxel);
      if (shift == kInfinity) break;

      totalShift += shift;
      if (minDistance <= totalShift) break;

      currentPoint += direction * (shift + shiftBonus);

      old = &candidates;
    }
    while (fVoxels.UpdateCurrentVoxel(currentPoint, direction, curVoxel));

    if (minCandidate < 0)
    {
      // No intersection found
      minDistance = 0;
      aConvex = false;
      Normal(aPoint, aNormalVector);
    }
    else
    {
      aConvex = (fExtremeFacets.find(fFacets[minCandidate])
              != fExtremeFacets.end());
    }
  }
  else
  {
    minDistance = DistanceToOutNoVoxels(aPoint, aDirection, aNormalVector,
                                        aConvex, aPstep);
  }
  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::
DistanceToInCandidates(const std::vector<G4int> &candidates,
                       const G4ThreeVector &aPoint,
                       const G4ThreeVector &direction) const
{
  G4int candidatesCount = candidates.size();
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal;

  G4double minDistance = kInfinity;   
  for (G4int i = 0 ; i < candidatesCount; ++i)
  {
    G4int candidate = candidates[i];
    G4VFacet &facet = *fFacets[candidate];
    if (facet.Intersect(aPoint,direction,false,dist,distFromSurface,normal))
    {
      //
      // Set minDist to the new distance to current facet if distFromSurface is
      // in positive direction and point is not at surface. If the point is
      // within 0.5*kCarTolerance of the surface, then force distance to be
      // zero and leave member function immediately (for efficiency), as
      // proposed by & credit to Akira Okumura.
      //
      if ( (distFromSurface > kCarToleranceHalf)
        && (dist >= 0.0) && (dist < minDistance))
      {
        minDistance  = dist;
      }
      else
      {
        if (-kCarToleranceHalf <= dist && dist <= kCarToleranceHalf)
        {
         return 0.0;
        }
        else if  (distFromSurface > -kCarToleranceHalf
               && distFromSurface <  kCarToleranceHalf)
        {
          minDistance = dist; 
        }
      }
    }
  }
  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double
G4TessellatedSolid::DistanceToInCore(const G4ThreeVector &aPoint,
                                     const G4ThreeVector &aDirection,
                                           G4double aPstep) const
{
  G4double minDistance;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    minDistance = kInfinity;
    G4ThreeVector currentPoint = aPoint;
    G4ThreeVector direction = aDirection.unit();
    G4double shift = fVoxels.DistanceToFirst(currentPoint, direction);
    if (shift == kInfinity) return shift;
    G4double shiftBonus = kCarTolerance;
    if (shift) 
      currentPoint += direction * (shift + shiftBonus);
    // if (!fVoxels.Contains(currentPoint))  return minDistance;
    G4double totalShift = shift;

    // G4SurfBits exclusion; // (1/*fVoxels.GetBitsPerSlice()*/);
    vector<G4int> curVoxel(3);

    fVoxels.GetVoxel(curVoxel, currentPoint);
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      const vector<G4int> &candidates = fVoxels.GetCandidates(curVoxel);
      if (candidates.size())
      {
        G4double distance=DistanceToInCandidates(candidates, aPoint, direction);
        if (minDistance > distance) minDistance = distance;
        if (distance < totalShift) break;
      }

      shift = fVoxels.DistanceToNext(currentPoint, direction, curVoxel);
      if (shift == kInfinity /*|| shift == 0*/) break;

      totalShift += shift;
      if (minDistance < totalShift) break;

      currentPoint += direction * (shift + shiftBonus);
    }
    while (fVoxels.UpdateCurrentVoxel(currentPoint, direction, curVoxel));
  }
  else
  {
    minDistance = DistanceToInNoVoxels(aPoint, aDirection, aPstep);
  }

  return minDistance;
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool
G4TessellatedSolid::CompareSortedVoxel(const std::pair<G4int, G4double> &l,
                                       const std::pair<G4int, G4double> &r)
{
  return l.second < r.second;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double
G4TessellatedSolid::MinDistanceFacet(const G4ThreeVector &p,
                                           G4bool simple,
                                           G4VFacet * &minFacet) const
{
  G4double minDist = kInfinity;

  G4int size = fVoxels.GetVoxelBoxesSize();
  vector<pair<G4int, G4double> > voxelsSorted(size);
  
  pair<G4int, G4double> info;

  for (G4int i = 0; i < size; ++i)
  {
    const G4VoxelBox &voxelBox = fVoxels.GetVoxelBox(i);

    G4ThreeVector pointShifted = p - voxelBox.pos;
    G4double safety = fVoxels.MinDistanceToBox(pointShifted, voxelBox.hlen);
    info.first = i;
    info.second = safety;
    voxelsSorted[i] = info;
  }

  std::sort(voxelsSorted.begin(), voxelsSorted.end(),
            &G4TessellatedSolid::CompareSortedVoxel);

  for (G4int i = 0; i < size; ++i)
  {
    const pair<G4int,G4double> &inf = voxelsSorted[i];
    G4double dist = inf.second;
    if (dist > minDist) break;

    const vector<G4int> &candidates = fVoxels.GetVoxelBoxCandidates(inf.first);
    G4int csize = candidates.size();
    for (G4int j = 0; j < csize; ++j)
    {
      G4int candidate = candidates[j];
      G4VFacet &facet = *fFacets[candidate];
      dist = simple ? facet.Distance(p,minDist)
                    : facet.Distance(p,minDist,false);
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
G4double G4TessellatedSolid::SafetyFromOutside (const G4ThreeVector &p,
                                                      G4bool aAccurate) const
{
#if G4SPECSDEBUG
  if ( Inside(p) == kInside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16) ;
    message << "Point p is already inside!?" << G4endl
      << "Position:"  << G4endl << G4endl
      << "p.x() = "   << p.x()/mm << " mm" << G4endl
      << "p.y() = "   << p.y()/mm << " mm" << G4endl
      << "p.z() = "   << p.z()/mm << " mm" << G4endl
      << "DistanceToOut(p) == " << DistanceToOut(p);
    message.precision(oldprc) ;
    G4Exception("G4TriangularFacet::DistanceToIn(p)",
                "GeomSolids1002", JustWarning, message);
  }
#endif

  G4double minDist;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    if (!aAccurate)
      return fVoxels.DistanceToBoundingBox(p);

    if (!OutsideOfExtent(p, kCarTolerance))
    {
      vector<G4int> startingVoxel(3);
      fVoxels.GetVoxel(startingVoxel, p);
      const vector<G4int> &candidates = fVoxels.GetCandidates(startingVoxel);
      if (candidates.size() == 0 && fInsides.GetNbits())
      {
        G4int index = fVoxels.GetPointIndex(p);
        if (fInsides[index]) return 0.;
      }
    }

    G4VFacet *facet;
    minDist = MinDistanceFacet(p, true, facet);
  }
  else
  {
    minDist = kInfinity;
    G4int size = fFacets.size();
    for (G4int i = 0; i < size; ++i)
    {
      G4VFacet &facet = *fFacets[i];
      G4double dist = facet.Distance(p,minDist);
      if (dist < minDist) minDist  = dist;
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double
G4TessellatedSolid::SafetyFromInside (const G4ThreeVector &p, G4bool) const
{  
#if G4SPECSDEBUG
  if ( Inside(p) == kOutside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16) ;
    message << "Point p is already outside!?" << G4endl
      << "Position:"  << G4endl << G4endl
      << "p.x() = "   << p.x()/mm << " mm" << G4endl
      << "p.y() = "   << p.y()/mm << " mm" << G4endl
      << "p.z() = "   << p.z()/mm << " mm" << G4endl
      << "DistanceToIn(p) == " << DistanceToIn(p);
    message.precision(oldprc) ;
    G4Exception("G4TriangularFacet::DistanceToOut(p)",
                "GeomSolids1002", JustWarning, message);
  }
#endif

  G4double minDist;

  if (OutsideOfExtent(p, kCarTolerance)) return 0.0;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    G4VFacet *facet;
    minDist = MinDistanceFacet(p, true, facet);
  }
  else
  {
    minDist = kInfinity;
    G4double dist = 0.0;
    G4int size = fFacets.size();
    for (G4int i = 0; i < size; ++i)
    {
      G4VFacet &facet = *fFacets[i];
      dist = facet.Distance(p,minDist);
      if (dist < minDist) minDist  = dist;
    }
  }
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4GeometryType GetEntityType() const;
//
// Provide identification of the class of an object
//
G4GeometryType G4TessellatedSolid::GetEntityType () const
{
  return fGeometryType;
}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream &G4TessellatedSolid::StreamInfo(std::ostream &os) const
{
  os << G4endl;
  os << "Geometry Type    = " << fGeometryType  << G4endl;
  os << "Number of facets = " << fFacets.size() << G4endl;

  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    os << "FACET #          = " << i + 1 << G4endl;
    G4VFacet &facet = *fFacets[i];
    facet.StreamInfo(os);
  }
  os << G4endl;

  return os;
}

///////////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4TessellatedSolid::Clone() const
{
  return new G4TessellatedSolid(*this);
}

///////////////////////////////////////////////////////////////////////////////
//
// EInside G4TessellatedSolid::Inside (const G4ThreeVector &p) const
//
// This method must return:
//    * kOutside if the point at offset p is outside the shape
//      boundaries plus kCarTolerance/2,
//    * kSurface if the point is <= kCarTolerance/2 from a surface, or
//    * kInside otherwise.
//
EInside G4TessellatedSolid::Inside (const G4ThreeVector &aPoint) const
{
  EInside location;

  if (fVoxels.GetCountOfVoxels() > 1)
  {
    location = InsideVoxels(aPoint);
  }
  else
  {
    location = InsideNoVoxels(aPoint);
  }
  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TessellatedSolid::SurfaceNormal(const G4ThreeVector& p) const
{
  G4ThreeVector n;
  Normal(p, n);
  return n;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an outside point p. The
// distance can be an underestimate.
//
G4double G4TessellatedSolid::DistanceToIn(const G4ThreeVector& p) const
{
  return SafetyFromOutside(p,false);
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToIn(const G4ThreeVector& p,
                                          const G4ThreeVector& v)const
{
  G4double dist = DistanceToInCore(p,v,kInfinity);
#ifdef G4SPECSDEBUG
  if (dist < kInfinity)
  {
    if (Inside(p + dist*v) != kSurface)
    {
      std::ostringstream message;
      message << "Invalid response from facet in solid '" << GetName() << "',"
              << G4endl
              << "at point: " << p <<  "and direction: " << v;
      G4Exception("G4TessellatedSolid::DistanceToIn(p,v)",
                  "GeomSolids1002", JustWarning, message);
    }
  }
#endif
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an inside
// point. The distance can be an underestimate.
//
G4double G4TessellatedSolid::DistanceToOut(const G4ThreeVector& p) const
{
  return SafetyFromInside(p,false);
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
//                        const G4bool calcNorm=false,
//                        G4bool *validNorm=0, G4ThreeVector *n=0);
//
// Return distance along the normalised vector v to the shape, from a
// point at an offset p inside or on the surface of the
// shape. Intersections with surfaces, when the point is not greater
// than kCarTolerance/2 from a surface, must be ignored.
//     If calcNorm is true, then it must also set validNorm to either
//     * true, if the solid lies entirely behind or on the exiting
//        surface. Then it must set n to the outwards normal vector
//        (the Magnitude of the vector is not defined).
//     * false, if the solid does not lie entirely behind or on the
//       exiting surface.
// If calcNorm is false, then validNorm and n are unused.
//
G4double G4TessellatedSolid::DistanceToOut(const G4ThreeVector& p,
                                           const G4ThreeVector& v,
                                           const G4bool calcNorm,
                                                 G4bool *validNorm,
                                                 G4ThreeVector *norm) const
{
  G4ThreeVector n;
  G4bool valid;

  G4double dist = DistanceToOutCore(p, v, n, valid);
  if (calcNorm)
  {
    *norm = n;
    *validNorm = valid;
  }
#ifdef G4SPECSDEBUG
  if (dist < kInfinity)
  {
    if (Inside(p + dist*v) != kSurface)
    {
      std::ostringstream message;
      message << "Invalid response from facet in solid '" << GetName() << "',"
              << G4endl
              << "at point: " << p <<  "and direction: " << v;
      G4Exception("G4TessellatedSolid::DistanceToOut(p,v,..)",
                  "GeomSolids1002", JustWarning, message);
    }
  }
#endif
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

///////////////////////////////////////////////////////////////////////////////
//
G4Polyhedron *G4TessellatedSolid::CreatePolyhedron () const
{
  G4int nVertices = fVertexList.size();
  G4int nFacets   = fFacets.size();
  G4PolyhedronArbitrary *polyhedron =
    new G4PolyhedronArbitrary (nVertices, nFacets);
  for (G4ThreeVectorList::const_iterator v= fVertexList.begin();
                                         v!=fVertexList.end(); ++v)
  {
    polyhedron->AddVertex(*v);
  }

  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    G4int v[4];
    G4int n = facet.GetNumberOfVertices();
    if (n > 4) n = 4;
    else if (n == 3) v[3] = 0;
    for (G4int j=0; j<n; ++j)
    {
      G4int k = facet.GetVertexIndex(j);
      v[j] = k+1;
    }
    polyhedron->AddFacet(v[0],v[1],v[2],v[3]);
  }
  polyhedron->SetReferences();  

  return (G4Polyhedron*) polyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron
//
G4Polyhedron* G4TessellatedSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
  {
    G4AutoLock l(&polyhedronMutex);
    delete fpPolyhedron;
    fpPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fpPolyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
// Get bounding box
//
void G4TessellatedSolid::BoundingLimits(G4ThreeVector& pMin,
                                        G4ThreeVector& pMax) const
{
  pMin = fMinExtent;
  pMax = fMaxExtent;

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4TessellatedSolid::BoundingLimits()",
                "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
//
G4bool
G4TessellatedSolid::CalculateExtent(const EAxis pAxis,
                                    const G4VoxelLimits& pVoxelLimit,
                                    const G4AffineTransform& pTransform,
                                          G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // The extent is calculated as cumulative extent of the pyramids
  // formed by facets and the center of the bounding box.
  //
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  G4ThreeVectorList base;
  G4ThreeVectorList apex(1);
  std::vector<const G4ThreeVectorList *> pyramid(2);
  pyramid[0] = &base;
  pyramid[1] = &apex;
  apex[0] = (bmin+bmax)*0.5;

  // main loop along facets
  pMin =  kInfinity;
  pMax = -kInfinity;
  for (G4int i=0; i<GetNumberOfFacets(); ++i)
  {
    G4VFacet* facet = GetFacet(i);
    if (std::abs((facet->GetSurfaceNormal()).dot(facet->GetVertex(0)-apex[0]))
        < kCarToleranceHalf) continue;

    G4int nv = facet->GetNumberOfVertices();
    base.resize(nv);
    for (G4int k=0; k<nv; ++k) { base[k] = facet->GetVertex(k); }

    G4double emin,emax;
    G4BoundingEnvelope benv(pyramid);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax)) continue;
    if (emin < pMin) pMin = emin;
    if (emax > pMax) pMax = emax;
    if (eminlim > pMin && emaxlim < pMax) break; // max possible extent
  }
  return (pMin < pMax);
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinXExtent () const
{
  return fMinExtent.x();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxXExtent () const
{
  return fMaxExtent.x();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinYExtent () const
{
  return fMinExtent.y();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxYExtent () const
{
  return fMaxExtent.y();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinZExtent () const
{
  return fMinExtent.z();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxZExtent () const
{
  return fMaxExtent.z();
}

///////////////////////////////////////////////////////////////////////////////
//
G4VisExtent G4TessellatedSolid::GetExtent () const
{
  return G4VisExtent (fMinExtent.x(), fMaxExtent.x(),
                      fMinExtent.y(), fMaxExtent.y(),
                      fMinExtent.z(), fMaxExtent.z());
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetCubicVolume ()
{
  if (fCubicVolume != 0.) return fCubicVolume;

  // For explanation of the following algorithm see:
  // https://en.wikipedia.org/wiki/Polyhedron#Volume
  // http://wwwf.imperial.ac.uk/~rn/centroid.pdf

  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    G4double area = facet.GetArea();
    G4ThreeVector unit_normal = facet.GetSurfaceNormal();
    fCubicVolume += area * (facet.GetVertex(0).dot(unit_normal));
  }
  fCubicVolume /= 3.;
  return fCubicVolume;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetSurfaceArea ()
{
  if (fSurfaceArea != 0.) return fSurfaceArea;

  G4int size = fFacets.size();
  for (G4int i = 0; i < size; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    fSurfaceArea += facet.GetArea();
  }
  return fSurfaceArea;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TessellatedSolid::GetPointOnSurface() const
{
  // Select randomly a facet and return a random point on it

  G4int i = (G4int) G4RandFlat::shoot(0., fFacets.size());
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
void G4TessellatedSolid::SetRandomVectors ()
{
  fRandir.resize(20);
  fRandir[0]  =
    G4ThreeVector(-0.9577428892113370, 0.2732676269591740, 0.0897405271949221);
  fRandir[1]  =
    G4ThreeVector(-0.8331264504940770,-0.5162067214954600,-0.1985722492445700);
  fRandir[2]  =
    G4ThreeVector(-0.1516671651108820, 0.9666292616127460, 0.2064580868390110);
  fRandir[3]  =
    G4ThreeVector( 0.6570250350323190,-0.6944539025883300, 0.2933460081893360);
  fRandir[4]  =
    G4ThreeVector(-0.4820456281280320,-0.6331060000098690,-0.6056474264406270);
  fRandir[5]  =
    G4ThreeVector( 0.7629032554236800 , 0.1016854697539910,-0.6384658864065180);
  fRandir[6]  =
    G4ThreeVector( 0.7689540409061150, 0.5034929891988220, 0.3939600142169160);
  fRandir[7]  =
    G4ThreeVector( 0.5765188359255740, 0.5997271636278330,-0.5549354566343150);
  fRandir[8]  =
    G4ThreeVector( 0.6660632777862070,-0.6362809868288380, 0.3892379937580790);
  fRandir[9]  =
    G4ThreeVector( 0.3824415020414780, 0.6541792713761380,-0.6525243125110690);
  fRandir[10] =
    G4ThreeVector(-0.5107726564526760, 0.6020905056811610, 0.6136760679616570);
  fRandir[11] =
    G4ThreeVector( 0.7459135439578050, 0.6618796061649330, 0.0743530220183488);
  fRandir[12] =
    G4ThreeVector( 0.1536405855311580, 0.8117477913978260,-0.5634359711967240);
  fRandir[13] =
    G4ThreeVector( 0.0744395301705579,-0.8707110101772920,-0.4861286795736560);
  fRandir[14] =
    G4ThreeVector(-0.1665874645185400, 0.6018553940549240,-0.7810369397872780);
  fRandir[15] =
    G4ThreeVector( 0.7766902003633100, 0.6014617505959970,-0.1870724331097450);
  fRandir[16] =
    G4ThreeVector(-0.8710128685847430,-0.1434320216603030,-0.4698551243971010);
  fRandir[17] =
    G4ThreeVector( 0.8901082092766820,-0.4388411398893870, 0.1229871120030100);
  fRandir[18] =
    G4ThreeVector(-0.6430417431544370,-0.3295938228697690, 0.6912779675984150);
  fRandir[19] =
    G4ThreeVector( 0.6331124368380410, 0.6306211461665000, 0.4488714875425340);

  fMaxTries = 20;
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::AllocatedMemoryWithoutVoxels()
{
  G4int base = sizeof(*this);
  base += fVertexList.capacity() * sizeof(G4ThreeVector);
  base += fRandir.capacity() * sizeof(G4ThreeVector);

  G4int limit = fFacets.size();
  for (G4int i = 0; i < limit; i++)
  {
    G4VFacet &facet = *fFacets[i];
    base += facet.AllocatedMemory();
  }

  std::set<G4VFacet *>::const_iterator beg, end, it;
  beg = fExtremeFacets.begin();
  end = fExtremeFacets.end();
  for (it = beg; it != end; it++)
  {
    G4VFacet &facet = *(*it);
    base += facet.AllocatedMemory();
  }
  return base;
}

///////////////////////////////////////////////////////////////////////////////
//
G4int G4TessellatedSolid::AllocatedMemory()
{
  G4int size = AllocatedMemoryWithoutVoxels();
  G4int sizeInsides = fInsides.GetNbytes();
  G4int sizeVoxels = fVoxels.AllocatedMemory();
  size += sizeInsides + sizeVoxels;
  return size;
}
