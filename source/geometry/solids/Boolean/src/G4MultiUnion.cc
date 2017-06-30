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
//
// $Id:$
//
// Implementation for G4MultiUnion class
//
// History:
//
// 06.04.17 G.Cosmo - Imported implementation in Geant4 for VecGeom migration
// 19.10.12 Marek Gayer - Original implementation from USolids module
// --------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "G4MultiUnion.hh"
#include "Randomize.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"
#include "G4AffineTransform.hh"
#include "G4DisplacedSolid.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

//______________________________________________________________________________
G4MultiUnion::G4MultiUnion(const G4String& name)
  : G4VSolid(name), fAccurate(false),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
  SetName(name);
  fSolids.clear();
  fTransformObjs.clear();
  fCubicVolume = 0;
  fSurfaceArea = 0;
  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
}

//______________________________________________________________________________
G4MultiUnion::~G4MultiUnion()
{
}

//______________________________________________________________________________
void G4MultiUnion::AddNode(G4VSolid& solid, G4Transform3D& trans)
{
  fSolids.push_back(&solid);
  fTransformObjs.push_back(trans);  // Store a local copy of transformations
}

//______________________________________________________________________________
G4VSolid* G4MultiUnion::Clone() const
{
  return new G4MultiUnion(*this);
}

// Copy constructor
//______________________________________________________________________________
G4MultiUnion::G4MultiUnion(const G4MultiUnion& rhs)
  : G4VSolid(rhs),fCubicVolume (rhs.fCubicVolume),
    fSurfaceArea (rhs.fSurfaceArea),
    kRadTolerance(rhs.kRadTolerance), fAccurate(false),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
}

// Fake default constructor for persistency
//______________________________________________________________________________
G4MultiUnion::G4MultiUnion( __void__& a )
  : G4VSolid(a), fCubicVolume (0.), fSurfaceArea (0.), kRadTolerance(0.),
    fAccurate(false), fRebuildPolyhedron(false), fpPolyhedron(0)
{
}

// Assignment operator
//______________________________________________________________________________
G4MultiUnion& G4MultiUnion::operator = (const G4MultiUnion& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  return *this;
}

//______________________________________________________________________________
G4double G4MultiUnion::GetCubicVolume()
{
  // Computes the cubic volume of the "G4MultiUnion" structure using
  // random points

  if (!fCubicVolume)
  {
    G4ThreeVector extentMin, extentMax, d, p, point;
    G4int inside = 0, generated;
    BoundingLimits(extentMin, extentMax);
    d = (extentMax - extentMin) / 2.;
    p = (extentMax + extentMin) / 2.;
    G4ThreeVector left = p - d;
    G4ThreeVector length = d * 2;
    for (generated = 0; generated < 10000; generated++)
    {
      G4ThreeVector rvec(G4UniformRand(), G4UniformRand(), G4UniformRand());
      point = left + G4ThreeVector(length.x()*rvec.x(),
                                   length.y()*rvec.y(),
                                   length.z()*rvec.z());
      if (Inside(point) != EInside::kOutside) inside++;
    }
    G4double vbox = (2 * d.x()) * (2 * d.y()) * (2 * d.z());
    fCubicVolume = inside * vbox / generated;
  }
  return fCubicVolume;
}

//______________________________________________________________________________
G4double
G4MultiUnion::DistanceToInNoVoxels(const G4ThreeVector& aPoint,
                                   const G4ThreeVector& aDirection) const
{
  G4ThreeVector direction = aDirection.unit();
  G4ThreeVector localPoint, localDirection;
  G4double minDistance = kInfinity;

  G4int numNodes = fSolids.size();
  for (G4int i = 0 ; i < numNodes ; ++i)
  {
    G4VSolid& solid = *fSolids[i];
    const G4Transform3D& transform = fTransformObjs[i];

    localPoint = GetLocalPoint(transform, aPoint);
    localDirection = GetLocalVector(transform, direction);

    G4double distance = solid.DistanceToIn(localPoint, localDirection);
    if (minDistance > distance) minDistance = distance;
  }
  return minDistance;
}

//______________________________________________________________________________
G4double G4MultiUnion::DistanceToInCandidates(const G4ThreeVector& aPoint,
                                              const G4ThreeVector& direction,
                                              std::vector<G4int>& candidates,
                                              G4SurfBits& bits) const
{
  G4int candidatesCount = candidates.size();
  G4ThreeVector localPoint, localDirection;

  G4double minDistance = kInfinity;
  for (G4int i = 0 ; i < candidatesCount; ++i)
  {
    G4int candidate = candidates[i];
    G4VSolid& solid = *fSolids[candidate];
    const G4Transform3D& transform = fTransformObjs[candidate];

    localPoint = GetLocalPoint(transform, aPoint);
    localDirection = GetLocalVector(transform, direction);
    G4double distance = solid.DistanceToIn(localPoint, localDirection);
    if (minDistance > distance) minDistance = distance;
    bits.SetBitNumber(candidate);
    if (minDistance == 0) break;
  }
  return minDistance;
}

// Algorithm note: we have to look also for all other objects in next voxels,
// if the distance is not shorter ... we have to do it because,
// for example for objects which starts in first voxel in which they
// do not collide with direction line, but in second it collides...
// The idea of crossing voxels would be still applicable,
// because this way we could exclude from the testing such solids,
// which were found that obviously are not good candidates, because
// they would return infinity
// But if distance is smaller than the shift to next voxel, we can return
// it immediately
//______________________________________________________________________________
G4double G4MultiUnion::DistanceToIn(const G4ThreeVector& aPoint,
                                    const G4ThreeVector& aDirection) const
{
  G4double minDistance = kInfinity;
  G4ThreeVector direction = aDirection.unit();
  G4double shift = fVoxels.DistanceToFirst(aPoint, direction);
  if (shift == kInfinity) return shift;

  G4ThreeVector currentPoint = aPoint;
  if (shift) currentPoint += direction * shift;

  G4SurfBits exclusion(fVoxels.GetBitsPerSlice());
  std::vector<G4int> candidates, curVoxel(3);
  fVoxels.GetVoxel(curVoxel, currentPoint);

  do
  {
    {
      if (fVoxels.GetCandidatesVoxelArray(curVoxel, candidates, &exclusion))
      {
        G4double distance = DistanceToInCandidates(aPoint, direction,
                                                   candidates, exclusion);
        if (minDistance > distance) minDistance = distance;
        if (distance < shift) break;
      }
    }
    shift = fVoxels.DistanceToNext(aPoint, direction, curVoxel);
  }
  while (minDistance > shift);
 
  return minDistance;
}

//______________________________________________________________________________
G4double G4MultiUnion::DistanceToOutNoVoxels(const G4ThreeVector& aPoint,
                                             const G4ThreeVector& aDirection,
                                             G4ThreeVector* aNormal) const
{
  // Computes distance from a point presumably outside the solid to the solid
  // surface. Ignores first surface if the point is actually inside.
  // Early return infinity in case the safety to any surface is found greater
  // than the proposed step aPstep.
  // The normal vector to the crossed surface is filled only in case the box
  // is crossed, otherwise aNormal->IsNull() is true.

  // algorithm:
  G4ThreeVector direction = aDirection.unit();
  G4ThreeVector localPoint, localDirection;
  G4int ignoredSolid = -1;
  G4double resultDistToOut = 0;
  G4ThreeVector currentPoint = aPoint;

  G4int numNodes = fSolids.size();
  for (G4int i = 0; i < numNodes; ++i)
  {
    if (i != ignoredSolid)
    {
      G4VSolid& solid = *fSolids[i];
      const G4Transform3D& transform = fTransformObjs[i];
      localPoint = GetLocalPoint(transform, currentPoint);
      localDirection = GetLocalVector(transform, direction);
      EInside location = solid.Inside(localPoint);
      if (location != EInside::kOutside)
      {
        G4double distance = solid.DistanceToOut(localPoint, localDirection,
                                                aNormal);
        if (distance < kInfinity)
        {
          if (resultDistToOut == kInfinity) resultDistToOut = 0;
          if (distance > 0)
          {
            currentPoint = GetGlobalPoint(transform, localPoint
                                          + distance*localDirection);
            resultDistToOut += distance;
            ignoredSolid = i; // skip the solid which we have just left
            i = -1; // force the loop to continue from 0
          }
        }
      }
    }
  }
  return resultDistToOut;
}

//______________________________________________________________________________
G4double G4MultiUnion::DistanceToOut(const G4ThreeVector& aPoint,
                                     const G4ThreeVector& aDirection,
                                     const G4bool /* calcNorm */,
                                     G4bool* /* validNorm */,
                                     G4ThreeVector* aNormal) const
{
  G4double distanceToOutVoxels = DistanceToOutVoxels(aPoint, aDirection,
                                                     aNormal);
  return distanceToOutVoxels;
}

//______________________________________________________________________________
G4double G4MultiUnion::DistanceToOutVoxels(const G4ThreeVector& aPoint,
                                           const G4ThreeVector& aDirection,
                                           G4ThreeVector* aNormal) const
{
  // Computes distance from a point presumably inside the solid to the solid
  // surface. Ignores first surface along each axis systematically (for points
  // inside or outside. Early returns zero in case the second surface is behind
  // the starting point.
  // o The proposed step is ignored.
  // o The normal vector to the crossed surface is always filled.

  // In the case the considered point is located inside the G4MultiUnion
  // structure, the treatments are as follows:
  //      - investigation of the candidates for the passed point
  //      - progressive moving of the point towards the surface, along the
  //        passed direction
  //      - processing of the normal

  G4ThreeVector direction = aDirection.unit();
  std::vector<G4int> candidates;
  G4double distance = 0;
  G4int numNodes = 2*fSolids.size();
  G4int count=0;

  if (fVoxels.GetCandidatesVoxelArray(aPoint, candidates))
  {
    // For normal case for which we presume the point is inside
    G4ThreeVector localPoint, localDirection, localNormal;
    G4ThreeVector currentPoint = aPoint;
    G4SurfBits exclusion(fVoxels.GetBitsPerSlice());
    G4bool notOutside;
    G4ThreeVector maxNormal;

    do
    {
      notOutside = false;

      G4double maxDistance = -kInfinity;
      G4int maxCandidate = 0;
      G4ThreeVector maxLocalPoint;

      G4int limit = candidates.size();
      for (G4int i = 0 ; i < limit ; ++i)
      {
        G4int candidate = candidates[i];
        // ignore the current component (that you just got out of) since
        // numerically the propagated point will be on its surface

        G4VSolid& solid = *fSolids[candidate];
        const G4Transform3D& transform = fTransformObjs[candidate];

        // The coordinates of the point are modified so as to fit the
        // intrinsic solid local frame:
        localPoint = GetLocalPoint(transform, currentPoint);

        // DistanceToOut at least for Trd sometimes return non-zero value
        // even from points that are outside. Therefore, this condition
        // must currently be here, otherwise it would not work.
        // But it means it would be slower.

        if (solid.Inside(localPoint) != EInside::kOutside)
        {
          notOutside = true;

          localDirection = GetLocalVector(transform, direction);

          // propagate with solid.DistanceToOut
          G4double shift = solid.DistanceToOut(localPoint, localDirection,
                                               false, 0, &localNormal);
          if (maxDistance < shift)
          {
            maxDistance = shift;
            maxCandidate = candidate;
            maxNormal = localNormal;
          }
        }
      }

      if (notOutside)
      {
        const G4Transform3D& transform = fTransformObjs[maxCandidate];

        // convert from local normal
        *aNormal = GetGlobalVector(transform, maxNormal);

        distance += maxDistance;
        currentPoint += maxDistance * direction;
        if(maxDistance == 0.) count++;

        // the current component will be ignored
        exclusion.SetBitNumber(maxCandidate);
        EInside location = InsideWithExclusion(currentPoint, &exclusion);

        // perform a Inside
        // it should be excluded current solid from checking
        // we have to collect the maximum distance from all given candidates.
        // such "maximum" candidate should be then used for finding next
        // candidates
        if (location == EInside::kOutside)
        {
          // else return cumulated distances to outside of the traversed
          // components
          break;
        }
        // if inside another component, redo 1 to 3 but add the next
        // DistanceToOut on top of the previous.

        // and fill the candidates for the corresponding voxel (just
        // exiting current component along direction)
        candidates.clear();

        fVoxels.GetCandidatesVoxelArray(currentPoint, candidates, &exclusion);
        exclusion.ResetBitNumber(maxCandidate);
      }
    }
    while  ((notOutside) && (count < numNodes));
  }

  return distance;
}

//______________________________________________________________________________
EInside G4MultiUnion::InsideWithExclusion(const G4ThreeVector& aPoint,
                                          G4SurfBits* exclusion) const
{
  // Classify point location with respect to solid:
  //  o eInside       - inside the solid
  //  o eSurface      - close to surface within tolerance
  //  o eOutside      - outside the solid

  // Hitherto, it is considered that only parallelepipedic nodes
  // can be added to the container

  // Implementation using voxelisation techniques:
  // ---------------------------------------------

  G4ThreeVector localPoint;
  EInside location = EInside::kOutside;

  std::vector<G4int> candidates;
  std::vector<G4MultiUnionSurface> surfaces;

  // TODO: test if it works well and if so measure performance
  // TODO: getPointIndex should not be used, instead GetVoxel + GetVoxelsIndex
  //       should be used
  // TODO: than pass result to GetVoxel further to GetCandidatesVoxelArray
  // TODO: eventually GetVoxel should be inlined here, early exit if any
  //       binary search is -1

  G4int limit = fVoxels.GetCandidatesVoxelArray(aPoint, candidates, exclusion);
  for (G4int i = 0 ; i < limit ; ++i)
  {
    G4int candidate = candidates[i];
    G4VSolid& solid = *fSolids[candidate];
    const G4Transform3D& transform = fTransformObjs[candidate];

    // The coordinates of the point are modified so as to fit the intrinsic
    // solid local frame:
    localPoint = GetLocalPoint(transform, aPoint);
    location = solid.Inside(localPoint);
    if (location == EInside::kInside) return EInside::kInside;
    else if (location == EInside::kSurface)
    {
      G4MultiUnionSurface surface;
      surface.point = localPoint;
      surface.solid = &solid;
      surfaces.push_back(surface);
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // Important comment: When two solids touch each other along a flat
  // surface, the surface points will be considered as kSurface, while points
  // located around will correspond to kInside (cf. G4UnionSolid)

  G4int size = surfaces.size();
  for (G4int i = 0; i < size - 1; ++i)
  {
    G4MultiUnionSurface& left = surfaces[i];
    for (G4int j = i + 1; j < size; ++j)
    {
      G4MultiUnionSurface& right = surfaces[j];
      G4ThreeVector n, n2;
      n = left.solid->SurfaceNormal(left.point);
      n2 = right.solid->SurfaceNormal(right.point);
      if ((n +  n2).mag2() < 1000 * kRadTolerance)
        return EInside::kInside;
    }
  }

  location = size ? EInside::kSurface : EInside::kOutside;

  return location;
}

//______________________________________________________________________________
EInside G4MultiUnion::Inside(const G4ThreeVector& aPoint) const
{
  // Classify point location with respect to solid:
  //  o eInside       - inside the solid
  //  o eSurface      - close to surface within tolerance
  //  o eOutside      - outside the solid

  // Hitherto, it is considered that only parallelepipedic nodes can be
  // added to the container

  // Implementation using voxelisation techniques:
  // ---------------------------------------------

  //  return InsideIterator(aPoint);

  EInside location = InsideWithExclusion(aPoint);
  return location;
}

//______________________________________________________________________________
EInside G4MultiUnion::InsideNoVoxels(const G4ThreeVector& aPoint) const
{
  G4ThreeVector localPoint;
  EInside location = EInside::kOutside;
  G4int countSurface = 0;

  G4int numNodes = fSolids.size();
  for (G4int i = 0 ; i < numNodes ; ++i)
  {
    G4VSolid& solid = *fSolids[i];
    G4Transform3D transform = GetTransformation(i);

    // The coordinates of the point are modified so as to fit the
    // intrinsic solid local frame:
    localPoint = GetLocalPoint(transform, aPoint);

    location = solid.Inside(localPoint);

    if (location == EInside::kSurface)
      countSurface++;

    if (location == EInside::kInside) return EInside::kInside;
  }
  if (countSurface != 0) return EInside::kSurface;
  return EInside::kOutside;
}

//______________________________________________________________________________
void G4MultiUnion::Extent(EAxis aAxis, G4double& aMin, G4double& aMax) const
{
  // Determines the bounding box for the considered instance of "UMultipleUnion"
  G4ThreeVector min, max;

  G4int numNodes = fSolids.size();
  for (G4int i = 0 ; i < numNodes ; ++i)
  {
    G4VSolid& solid = *fSolids[i];
    G4Transform3D transform = GetTransformation(i);
    solid.BoundingLimits(min, max);
   
    TransformLimits(min, max, transform);
    
    if (i == 0)
    {
      switch (aAxis)
      {
        case kXAxis:
          aMin = min.x();
          aMax = max.x();
          break;
        case kYAxis:
          aMin = min.y();
          aMax = max.y();
          break;
        case kZAxis:
          aMin = min.z();
          aMax = max.z();
          break;
        default:
          break;
      }
    }
    else
    {
      // Determine the min/max on the considered axis:
      switch (aAxis)
      {
        case kXAxis:
          if (min.x() < aMin)
            aMin = min.x();
          if (max.x() > aMax)
            aMax = max.x();
          break;
        case kYAxis:
          if (min.y() < aMin)
            aMin = min.y();
          if (max.y() > aMax)
            aMax = max.y();
          break;
        case kZAxis:
          if (min.z() < aMin)
            aMin = min.z();
          if (max.z() > aMax)
            aMax = max.z();
          break;
        default:
          break;
      }
    }
  }
}

//______________________________________________________________________________
void G4MultiUnion::BoundingLimits(G4ThreeVector& aMin,
                                  G4ThreeVector& aMax) const
{
   Extent(kXAxis, aMin[0], aMax[0]);
   Extent(kYAxis, aMin[1], aMax[1]);
   Extent(kZAxis, aMin[2], aMax[2]);
}

//______________________________________________________________________________
G4bool
G4MultiUnion::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
                                    G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

//______________________________________________________________________________
G4ThreeVector G4MultiUnion::SurfaceNormal(const G4ThreeVector& aPoint) const
{
  // Computes the localNormal on a surface and returns it as a unit vector.
  // Must return a valid vector. (even if the point is not on the surface).
  //
  //   On an edge or corner, provide an average localNormal of all facets within
  //   tolerance
  //   NOTE: the tolerance value used in here is not yet the global surface
  //         tolerance - we will have to revise this value - TODO

  std::vector<G4int> candidates;
  G4ThreeVector localPoint, normal, localNormal;
  G4double safety = kInfinity;
  G4int node = 0;

  ///////////////////////////////////////////////////////////////////////////
  // Important comment: Cases for which the point is located on an edge or
  // on a vertice remain to be treated

  // determine weather we are in voxel area
  if (fVoxels.GetCandidatesVoxelArray(aPoint, candidates))
  {
    G4int limit = candidates.size();
    for (G4int i = 0 ; i < limit ; ++i)
    {
      G4int candidate = candidates[i];
      const G4Transform3D& transform = fTransformObjs[candidate];

      // The coordinates of the point are modified so as to fit the intrinsic
      // solid local frame:
      localPoint = GetLocalPoint(transform, aPoint);
      G4VSolid& solid = *fSolids[candidate];
      EInside location = solid.Inside(localPoint);

      if (location == EInside::kSurface)
      {
        // normal case when point is on surface, we pick first solid
        normal = GetGlobalVector(transform, solid.SurfaceNormal(localPoint));
        return normal.unit();
      }
      else
      {
        // collect the smallest safety and remember solid node
        G4double s = (location == EInside::kInside)
                   ? solid.DistanceToOut(localPoint)
                   : solid.DistanceToIn(localPoint);
        if (s < safety)
        {
          safety = s;
          node = candidate;
        }
      }
    }
    // on none of the solids, the point was not on the surface
    G4VSolid& solid = *fSolids[node];
    const G4Transform3D& transform = fTransformObjs[node];
    localPoint = GetLocalPoint(transform, aPoint);

    normal = GetGlobalVector(transform, solid.SurfaceNormal(localPoint));
    return normal.unit();
  }
  else
  {
    // for the case when point is certainly outside:

    // find a solid in union with the smallest safety
    node = SafetyFromOutsideNumberNode(aPoint, safety);
    G4VSolid& solid = *fSolids[node];

    const G4Transform3D& transform = fTransformObjs[node];
    localPoint = GetLocalPoint(transform, aPoint);

    // evaluate normal for point at this found solid
    // and transform multi-union coordinates
    normal = GetGlobalVector(transform, solid.SurfaceNormal(localPoint));

    return normal.unit();
  }
}

//______________________________________________________________________________
G4double G4MultiUnion::DistanceToOut(const G4ThreeVector& point) const
{
  // Estimates isotropic distance to the surface of the solid. This must
  // be either accurate or an underestimate.
  // Two modes: - default/fast mode, sacrificing accuracy for speed
  //            - "precise" mode,  requests accurate value if available.

  std::vector<G4int> candidates;
  G4ThreeVector localPoint;
  G4double safetyMin = kInfinity;

  // In general, the value return by DistanceToIn(p) will not be the exact
  // but only an undervalue (cf. overlaps)
  fVoxels.GetCandidatesVoxelArray(point, candidates);

  G4int limit = candidates.size();
  for (G4int i = 0; i < limit; ++i)
  {
    G4int candidate = candidates[i];

    // The coordinates of the point are modified so as to fit the intrinsic
    // solid local frame:
    const G4Transform3D& transform = fTransformObjs[candidate];
    localPoint = GetLocalPoint(transform, point);
    G4VSolid& solid = *fSolids[candidate];
    if (solid.Inside(localPoint) == EInside::kInside)
    {
      G4double safety = solid.DistanceToOut(localPoint);
      if (safetyMin > safety) safetyMin = safety;
    }
  }
  if (safetyMin == kInfinity) safetyMin = 0; // we are not inside

  return safetyMin;
}

//______________________________________________________________________________
G4double G4MultiUnion::DistanceToIn(const G4ThreeVector& point) const
{
  // Estimates the isotropic safety from a point outside the current solid to
  // any of its surfaces. The algorithm may be accurate or should provide a fast
  // underestimate.

  if (!fAccurate)  { return fVoxels.DistanceToBoundingBox(point); }

  const std::vector<G4VoxelBox>& boxes = fVoxels.GetBoxes();
  G4double safetyMin = kInfinity;
  G4ThreeVector localPoint;

  G4int numNodes = fSolids.size();
  for (G4int j = 0; j < numNodes; ++j)
  {
    G4ThreeVector dxyz;
    if (j > 0)
    {
      const G4ThreeVector& pos = boxes[j].pos;
      const G4ThreeVector& hlen = boxes[j].hlen;
      for (G4int i = 0; i <= 2; ++i)
        // distance to middle point - hlength => distance from point to border
        // of x,y,z
        if ((dxyz[i] = std::abs(point[i] - pos[i]) - hlen[i]) > safetyMin)
          continue;

      G4double d2xyz = 0.;
      for (G4int i = 0; i <= 2; ++i)
        if (dxyz[i] > 0) d2xyz += dxyz[i] * dxyz[i];

      // minimal distance is at least this, but could be even higher. therefore,
      // we can stop if previous was already lower, let us check if it does any
      // chance to be better tha previous values...
      if (d2xyz >= safetyMin * safetyMin)
      {
        continue;
      }
    }
    const G4Transform3D& transform = fTransformObjs[j];
    localPoint = GetLocalPoint(transform, point);
    G4VSolid& solid = *fSolids[j];

    G4double safety = solid.DistanceToIn(localPoint);
    if (safety <= 0) return safety;
      // it was detected, that the point is not located outside
    if (safetyMin > safety) safetyMin = safety;
  }
  return safetyMin;
}

//______________________________________________________________________________
G4double G4MultiUnion::GetSurfaceArea()
{
  if (!fSurfaceArea)
  {
    fSurfaceArea = EstimateSurfaceArea(1000000, 0.001);
  }
  return fSurfaceArea;
}

//______________________________________________________________________________
void G4MultiUnion::Voxelize()
{
  fVoxels.Voxelize(fSolids, fTransformObjs);
}

//______________________________________________________________________________
G4int G4MultiUnion::SafetyFromOutsideNumberNode(const G4ThreeVector& aPoint,
                                                G4double& safetyMin) const
{
  // Method returning the closest node from a point located outside a
  // G4MultiUnion.
  // This is used to compute the normal in the case no candidate has been found.

  const std::vector<G4VoxelBox>& boxes = fVoxels.GetBoxes();
  safetyMin = kInfinity;
  G4int safetyNode = 0;
  G4ThreeVector localPoint;

  G4int numNodes = fSolids.size();
  for (G4int i = 0; i < numNodes; ++i)
  {
    G4double d2xyz = 0.;
    G4double dxyz0 = std::abs(aPoint.x() - boxes[i].pos.x()) - boxes[i].hlen.x();
    if (dxyz0 > safetyMin) continue;
    G4double dxyz1 = std::abs(aPoint.y() - boxes[i].pos.y()) - boxes[i].hlen.y();
    if (dxyz1 > safetyMin) continue;
    G4double dxyz2 = std::abs(aPoint.z() - boxes[i].pos.z()) - boxes[i].hlen.z();
    if (dxyz2 > safetyMin) continue;

    if (dxyz0 > 0) d2xyz += dxyz0 * dxyz0;
    if (dxyz1 > 0) d2xyz += dxyz1 * dxyz1;
    if (dxyz2 > 0) d2xyz += dxyz2 * dxyz2;
    if (d2xyz >= safetyMin * safetyMin) continue;

    G4VSolid& solid = *fSolids[i];
    const G4Transform3D& transform = fTransformObjs[i];
    localPoint = GetLocalPoint(transform, aPoint);
    fAccurate = true;
    G4double safety = solid.DistanceToIn(localPoint);
    fAccurate = false;
    if (safetyMin > safety)
    {
      safetyMin = safety;
      safetyNode = i;
    }
  }
  return safetyNode;
}

//______________________________________________________________________________
void G4MultiUnion::TransformLimits(G4ThreeVector& min, G4ThreeVector& max,
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

// Stream object contents to an output stream
//______________________________________________________________________________
std::ostream& G4MultiUnion::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "                *** Dump for solid - " << GetName() << " ***\n"
     << "                ===================================================\n"
     << " Solid type: G4MultiUnion\n"
     << " Parameters: \n";
     G4int numNodes = fSolids.size();
     for (G4int i = 0 ; i < numNodes ; ++i)
     {
      G4VSolid& solid = *fSolids[i];
      solid.StreamInfo(os);
      const G4Transform3D& transform = fTransformObjs[i];
      os << " Translation is " << transform.getTranslation() << " \n";
      os << " Rotation is :" << " \n";
      os << " " << transform.getRotation() << "\n";
     }
  os << "             \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//______________________________________________________________________________
G4ThreeVector G4MultiUnion::GetPointOnSurface() const
{
  G4ThreeVector point;

  G4long size = fSolids.size();

  do
  {
    G4long rnd = G4RandFlat::shootInt(G4long(0), size);
    G4VSolid& solid = *fSolids[rnd];
    point = solid.GetPointOnSurface();
    const G4Transform3D& transform = fTransformObjs[rnd];
    point = GetGlobalPoint(transform, point);
  }
  while (Inside(point) != EInside::kSurface);

  return point;
}

//______________________________________________________________________________
void 
G4MultiUnion::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

//______________________________________________________________________________
G4Polyhedron* G4MultiUnion::CreatePolyhedron() const
{
  HepPolyhedronProcessor processor;
  HepPolyhedronProcessor::Operation operation = HepPolyhedronProcessor::UNION;

  G4VSolid* solidA = GetSolid(0);
  const G4Transform3D transform0=GetTransformation(0);
  G4DisplacedSolid dispSolidA("placedA",solidA,transform0);

  G4Polyhedron* top = new G4Polyhedron(*dispSolidA.GetPolyhedron());
    
  for(G4int i=1; i<GetNumberOfSolids(); ++i)
  {
    G4VSolid* solidB = GetSolid(i);
    const G4Transform3D transform=GetTransformation(i);
    G4DisplacedSolid dispSolidB("placedB",solidB,transform);
    G4Polyhedron* operand = dispSolidB.GetPolyhedron();
    processor.push_back (operation, *operand);
  }
   
  if (processor.execute(*top)) { return top; }
  else { return 0; } 
}

//______________________________________________________________________________
G4Polyhedron* G4MultiUnion::GetPolyhedron() const
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
