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
// UMultiUnion
//
// 19.10.12 Marek Gayer
// --------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "UVoxelizer.hh"
#include "UMultiUnion.hh"
#include "UUtils.hh"

#include "UBox.hh"


using namespace std;

//______________________________________________________________________________
UMultiUnion::UMultiUnion(const std::string& name)
{
  SetName(name);
  fSolids.clear();
  fTransformObjs.clear();
  fCubicVolume = 0;
  fSurfaceArea = 0;
}

//______________________________________________________________________________
UMultiUnion::~UMultiUnion()
{
}

//______________________________________________________________________________
void UMultiUnion::AddNode(VUSolid& solid, UTransform3D& trans)
{
  fSolids.push_back(&solid);
  fTransformObjs.push_back(trans);  // Store a local copy of transformations
}

//______________________________________________________________________________
VUSolid* UMultiUnion:: Clone() const
{
  return new UMultiUnion(*this);
}
//
// Copy constructor
//______________________________________________________________________________
UMultiUnion::UMultiUnion(const UMultiUnion& rhs)
  : VUSolid(rhs),fCubicVolume (rhs.fCubicVolume),
    fSurfaceArea (rhs.fSurfaceArea)
{
}

//
// Assignment operator
//______________________________________________________________________________
UMultiUnion& UMultiUnion::operator = (const UMultiUnion& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  VUSolid::operator=(rhs);

  return *this;
}

//______________________________________________________________________________
double UMultiUnion::Capacity()
{
  // Capacity computes the cubic volume of the "UMultiUnion" structure using random points

  // Random initialization:
  //  srand((unsigned int)time(NULL));
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    UVector3 extentMin, extentMax, d, p, point;
    int inside = 0, generated;
    Extent(extentMin, extentMax);
    d = (extentMax - extentMin) / 2.;
    p = (extentMax + extentMin) / 2.;
    UVector3 left = p - d;
    UVector3 length = d * 2;
    for (generated = 0; generated < 10000; generated++)
    {
      UVector3 random(rand(), rand(), rand());
      point = left + length.MultiplyByComponents(random / RAND_MAX);
      if (Inside(point) != eOutside) inside++;
    }
    double vbox = (2 * d.x()) * (2 * d.y()) * (2 * d.z());
    fCubicVolume = inside * vbox / generated;
  }
  return fCubicVolume;
}

//______________________________________________________________________________
void UMultiUnion::ComputeBBox(UBBox* /*aBox*/, bool /*aStore*/)
{
  // Computes bounding box.
  cout << "ComputeBBox - Not implemented" << endl;
}

//______________________________________________________________________________
double UMultiUnion::DistanceToInNoVoxels(const UVector3& aPoint, const UVector3& aDirection, // UVector3 &aNormal,
                                         double aPstep) const
{
  UVector3 direction = aDirection.Unit();
  UVector3 localPoint, localDirection;
  double minDistance = UUtils::kInfinity;

  int numNodes = fSolids.size();
  for (int i = 0 ; i < numNodes ; ++i)
  {
    VUSolid& solid = *fSolids[i];
    const UTransform3D& transform = fTransformObjs[i];

    localPoint = transform.LocalPoint(aPoint);
    localDirection = transform.LocalVector(direction);

    double distance = solid.DistanceToIn(localPoint, localDirection, aPstep);
    if (minDistance > distance) minDistance = distance;
  }
  return minDistance;
}

//______________________________________________________________________________
double UMultiUnion::DistanceToInCandidates(const UVector3& aPoint, const UVector3& direction, double aPstep, std::vector<int>& candidates, UBits& bits) const
{
  int candidatesCount = candidates.size();
  UVector3 localPoint, localDirection;

  double minDistance = UUtils::kInfinity;
  for (int i = 0 ; i < candidatesCount; ++i)
  {
    int candidate = candidates[i];
    VUSolid& solid = *fSolids[candidate];
    const UTransform3D& transform = fTransformObjs[candidate];

    localPoint = transform.LocalPoint(aPoint);
    localDirection = transform.LocalVector(direction);
    double distance = solid.DistanceToIn(localPoint, localDirection, aPstep);
    if (minDistance > distance) minDistance = distance;
    bits.SetBitNumber(candidate);
    if (minDistance == 0) break;
    
  }
  return minDistance;
}

// Algorithm note: we have to look also for all other objects in next voxels, if the distance is not shorter ... we have to do it because,
// for example for objects which starts in first voxel in which they
// do not collide with direction line, but in second it collides...
// The idea of crossing voxels would be still applicable,
// because this way we could exclude from the testing such solids,
// which were found that obviously are not good candidates, because
// they would return infinity
// But if distance is smaller than the shift to next voxel, we can return it immediately
//______________________________________________________________________________
double UMultiUnion::DistanceToIn(const UVector3& aPoint,
                                 const UVector3& aDirection, double aPstep) const
{
  //return DistanceToInNoVoxels(aPoint, aDirection, aPstep);

#ifdef DEBUG
  double distanceToInNoVoxels = DistanceToInNoVoxels(aPoint, aDirection, aPstep);
#endif

  double minDistance = UUtils::kInfinity;
  UVector3 direction = aDirection.Unit();
  double shift = fVoxels.DistanceToFirst(aPoint, direction);
  if (shift == UUtils::kInfinity) return shift;

  UVector3 currentPoint = aPoint;
  if (shift) currentPoint += direction * shift;
  //    if (!fVoxels.Contains(currentPoint))
  //      return minDistance;

  UBits exclusion(fVoxels.GetBitsPerSlice());
  vector<int> candidates, curVoxel(3);
  fVoxels.GetVoxel(curVoxel, currentPoint);

  do
  {
//    if (!fVoxels.Empty().GetNbits() || !fVoxels.IsEmpty(fVoxels.GetVoxelsIndex(curVoxel)))
    {
      if (fVoxels.GetCandidatesVoxelArray(curVoxel, candidates, &exclusion))
      {
        double distance = DistanceToInCandidates(aPoint, direction, aPstep, candidates, exclusion);
        if (minDistance > distance) minDistance = distance;
        if (distance < shift) break;
      }
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
 
  return minDistance;
}

//______________________________________________________________________________
double UMultiUnion::DistanceToOutNoVoxels(const UVector3& aPoint, const UVector3& aDirection,
                                          UVector3& aNormal,
                                          bool&     convex,
                                          double   /*aPstep*/) const
{
  // Computes distance from a point presumably outside the solid to the solid
  // surface. Ignores first surface if the point is actually inside. Early return
  // infinity in case the safety to any surface is found greater than the proposed
  // step aPstep.
  // The normal vector to the crossed surface is filled only in case the box is
  // crossed, otherwise aNormal.IsNull() is true.

  // algorithm:
  UVector3 direction = aDirection.Unit();
  UVector3 localPoint, localDirection;
  int ignoredSolid = -1;
  double resultDistToOut = 0; // UUtils::kInfinity;
  UVector3 currentPoint = aPoint;

  int numNodes = fSolids.size();
  for (int i = 0; i < numNodes; ++i)
  {
    if (i != ignoredSolid)
    {
      VUSolid& solid = *fSolids[i];
      const UTransform3D& transform = fTransformObjs[i];
      localPoint = transform.LocalPoint(currentPoint);
      localDirection = transform.LocalVector(direction);
      VUSolid::EnumInside location = solid.Inside(localPoint);
      if (location != eOutside)
      {
        double distance = solid.DistanceToOut(localPoint, localDirection, aNormal, convex);
        if (distance < UUtils::kInfinity)
        {
          if (resultDistToOut == UUtils::kInfinity) resultDistToOut = 0;
          if (distance > 0)
          {
            currentPoint = transform.GlobalPoint(localPoint + distance * localDirection);
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
double UMultiUnion::DistanceToOut(const UVector3& aPoint, const UVector3& aDirection,
                                  UVector3& aNormal,
                                  bool&     convex,
                                  double   aPstep) const
{
  //return DistanceToOutNoVoxels(aPoint, aDirection, aNormal, convex, aPstep);

#ifdef DEBUG
  double distanceToOutNoVoxels = DistanceToOutNoVoxels(aPoint, aDirection, aNormal, convex, aPstep);
#endif

  double distanceToOutVoxels = DistanceToOutVoxels(aPoint, aDirection, aNormal, convex, aPstep);
 
#ifdef DEBUG
  if (std::abs(distanceToOutVoxels - distanceToOutNoVoxels) > VUSolid::Tolerance())
  {
    // distanceToOutVoxels = distanceToOutVoxels;
    Inside(aPoint);
  }
  //    return distanceToOutNoVoxels;
#endif
  return distanceToOutVoxels;
}

//______________________________________________________________________________
double UMultiUnion::DistanceToOutVoxels(const UVector3& aPoint, const UVector3& aDirection,
                                        UVector3& aNormal,
                                        bool&     /*convex*/,
                                        double   /*aPstep*/) const
{
  // Computes distance from a point presumably inside the solid to the solid
  // surface. Ignores first surface along each axis systematically (for points
  // inside or outside. Early returns zero in case the second surface is behind
  // the starting point.
  // o The proposed step is ignored.
  // o The normal vector to the crossed surface is always filled.

  // In the case the considered point is located inside the UMultiUnion structure,
  // the treatments are as follows:
  //      - investigation of the candidates for the passed point
  //      - progressive moving of the point towards the surface, along the passed direction
  //      - processing of the normal

  UVector3 direction = aDirection.Unit();
  vector<int> candidates;
  double distance = 0;
  int numNodes = 2*fSolids.size();
  int count=0;

  if (fVoxels.GetCandidatesVoxelArray(aPoint, candidates))
  {
    // For normal case for which we presume the point is inside
    UVector3 localPoint, localDirection, localNormal;
    UVector3 currentPoint = aPoint;
    UBits exclusion(fVoxels.GetBitsPerSlice());
    bool notOutside;
    UVector3 maxNormal;

    do
    {
      notOutside = false;

      double maxDistance = -UUtils::kInfinity;
      int maxCandidate = 0;
      UVector3 maxLocalPoint;

      int limit = candidates.size();
      for (int i = 0 ; i < limit ; ++i)
      {
        int candidate = candidates[i];
        // ignore the current component (that you just got out of) since numerically the propagated point will be on its surface

        VUSolid& solid = *fSolids[candidate];
        const UTransform3D& transform = fTransformObjs[candidate];

        // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
        localPoint = transform.LocalPoint(currentPoint);

        // DistanceToOut at least for Trd sometimes return non-zero value even from points that are outside. Therefore, this condition must currently be here, otherwise it would not work. But it means it would be slower.

        if (solid.Inside(localPoint) != eOutside)
        {
          notOutside = true;

          localDirection = transform.LocalVector(direction);
          // propagate with solid.DistanceToOut
          bool convex;

          double shift = solid.DistanceToOut(localPoint, localDirection, localNormal, convex);

//          if (shift != 0)
//          {
//            if(solid.Inside(localPoint) == eOutside)
//              shift = shift;
//            notOutside = true;

          if (maxDistance < shift)
          {
            maxDistance = shift;
            maxCandidate = candidate;
            maxNormal = localNormal;
//              maxLocalDirection = localDirection;
          }
//          }
        }
      }

      if (notOutside)
      {
        const UTransform3D& transform = fTransformObjs[maxCandidate];
//        localPoint = transform.LocalPoint(currentPoint);
        // convert from local normal
        aNormal = transform.GlobalVector(maxNormal);

        distance += maxDistance;
//        currentPoint = transform.GlobalPoint(localPoint+maxDistance*maxLocalDirection);
        currentPoint += maxDistance * direction;
         if(maxDistance == 0.)count++;
        // the current component will be ignored
        exclusion.SetBitNumber(maxCandidate);
        VUSolid::EnumInside location = InsideWithExclusion(currentPoint, &exclusion);

        // perform a Inside
        // it should be excluded current solid from checking
        // we have to collect the maximum distance from all given candidates. such "maximum" candidate should be then used for finding next candidates
        if (location == eOutside)
        {
          // else return cumulated distances to outside of the traversed components
          break;
        }
        // if inside another component, redo 1 to 3 but add the next DistanceToOut on top of the previous.

        // and fill the candidates for the corresponding voxel (just exiting current component along direction)
        candidates.clear();

        fVoxels.GetCandidatesVoxelArray(currentPoint, candidates, &exclusion);
        exclusion.ResetBitNumber(maxCandidate);
      }
    }
    while  ((notOutside) && (count < numNodes));

//    if (distance != -UUtils::kInfinity)
//      return distance;
  }

  return distance;
}


struct USurface
{
  UVector3 point;
  VUSolid* solid;
};

//______________________________________________________________________________
VUSolid::EnumInside UMultiUnion::InsideWithExclusion(const UVector3& aPoint, UBits* exclusion) const
{
  // Classify point location with respect to solid:
  //  o eInside       - inside the solid
  //  o eSurface      - close to surface within tolerance
  //  o eOutside      - outside the solid

  // Hitherto, it is considered that:
  //        - only parallelepipedic nodes can be added to the container

  // Implementation using voxelisation techniques:
  // ---------------------------------------------

  UVector3 localPoint;
  VUSolid::EnumInside location = eOutside;

  vector<int> candidates;
  vector<USurface> surfaces;

  // TODO: test if it works well and if so measure performance
  // TODO: getPointIndex should not be used, instead GetVoxel + GetVoxelsIndex should be used
  // TODO: than pass result to GetVoxel further to GetCandidatesVoxelArray
  // TODO: eventually GetVoxel should be inlined here, early exit if any binary search is -1

  // the code bellow makes it currently only slower
  //  if (!fVoxels.empty.GetNbits() || !fVoxels.empty[fVoxels.GetPointIndex(aPoint)])
  {
    int limit = fVoxels.GetCandidatesVoxelArray(aPoint, candidates, exclusion);
    for (int i = 0 ; i < limit ; ++i)
    {
      int candidate = candidates[i];
      VUSolid& solid = *fSolids[candidate];
      const UTransform3D& transform = fTransformObjs[candidate];

      // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
      localPoint = transform.LocalPoint(aPoint);
      location = solid.Inside(localPoint);
      if (location == eInside) return eInside;
      else if (location == eSurface)
      {
        USurface surface;
        surface.point = localPoint;
        surface.solid = &solid;
        surfaces.push_back(surface);
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////
  // Important comment: When two solids touch each other along a flat
  // surface, the surface points will be considered as eSurface, while points
  // located around will correspond to eInside (cf. G4UnionSolid in GEANT4)

  int size = surfaces.size();
  for (int i = 0; i < size - 1; ++i)
  {
    USurface& left = surfaces[i];
    for (int j = i + 1; j < size; ++j)
    {
      USurface& right = surfaces[j];
      UVector3 n, n2;
      left.solid->Normal(left.point, n);
      right.solid->Normal(right.point, n2);
      if ((n +  n2).Mag2() < 1000 * frTolerance)
        return eInside;
    }
  }

  location = size ? eSurface : eOutside;

  return location;
}

/*
//______________________________________________________________________________
VUSolid::EnumInside UMultiUnion::InsideWithExclusion(const UVector3 &aPoint, UBits *exclusion) const
{
  // Classify point location with respect to solid:
  //  o eInside       - inside the solid
  //  o eSurface      - close to surface within tolerance
  //  o eOutside      - outside the solid

  // Hitherto, it is considered that:
  //        - only parallelepipedic nodes can be added to the container

  // Implementation using voxelisation techniques:
  // ---------------------------------------------

  UVector3 localPoint;
  VUSolid::EnumInside location = eOutside;
  bool surface = false;

  // TODO: test if it works well and if so measure performance
  // TODO: getPointIndex should not be used, instead GetVoxel + GetVoxelsIndex should be used
  // TODO: than pass result to GetVoxel further to GetCandidatesVoxelArray
  // TODO: eventually GetVoxel should be inlined here, early exit if any binary search is -1

  // the code bellow makes it currently only slower
  //  if (!fVoxels.empty.GetNbits() || !fVoxels.empty[fVoxels.GetPointIndex(aPoint)])
  {
    vector<int> curVoxel(3);
    if (fVoxels.GetPointVoxel(aPoint, curVoxel))
    {
      int index = fVoxels.GetVoxelsIndex(curVoxel);
      if (!fVoxels.Empty()[index])
      {
    //    int limit = fVoxels.GetCandidatesVoxelArray(aPoint, candidates, exclusion);
        const vector<int> &candidates = fVoxels.GetCandidates(curVoxel);
        int limit = candidates.size();
        for(int i = 0 ; i < limit ; ++i)
        {
          int candidate = candidates[i];
          VUSolid &solid = *solids[candidate];
          UTransform3D &transform = *transforms[candidate];

          // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
          localPoint = transform.LocalPoint(aPoint);
          location = solid.Inside(localPoint);
          if(location == eSurface) surface = true;

          if(location == eInside) return eInside;
        }
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////////
  // Important comment: When two solids touch each other along a flat
  // surface, the surface points will be considered as eSurface, while points
  // located around will correspond to eInside (cf. G4UnionSolid in GEANT4)
  location = surface ? eSurface : eOutside;

  return location;
}
*/

#define DEBU

//______________________________________________________________________________
VUSolid::EnumInside UMultiUnion::Inside(const UVector3& aPoint) const
{
  // Classify point location with respect to solid:
  //  o eInside       - inside the solid
  //  o eSurface      - close to surface within tolerance
  //  o eOutside      - outside the solid

  // Hitherto, it is considered that:
  //        - only parallelepipedic nodes can be added to the container

  // Implementation using voxelisation techniques:
  // ---------------------------------------------

  //  return InsideIterator(aPoint);

#ifdef DEBUG
  VUSolid::EnumInside insideNoVoxels = InsideNoVoxels(aPoint);
#endif

  VUSolid::EnumInside location = InsideWithExclusion(aPoint);

#ifdef DEBUG
  if (location != insideNoVoxels)
    location = insideNoVoxels; // you can place a breakpoint here
#endif
  return location;
}

//______________________________________________________________________________
VUSolid::EnumInside UMultiUnion::InsideNoVoxels(const UVector3& aPoint) const
{
  UVector3 localPoint;
  VUSolid::EnumInside location = eOutside;
  int countSurface = 0;

  int numNodes = fSolids.size();
  for (int i = 0 ; i < numNodes ; ++i)
  {
    VUSolid& solid = *fSolids[i];
    UTransform3D transform = GetTransformation(i);

    // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
    localPoint = transform.LocalPoint(aPoint);

    location = solid.Inside(localPoint);

    if (location == eSurface)
      countSurface++;

    if (location == eInside) return eInside;
  }
  if (countSurface != 0) return eSurface;
  return eOutside;
}

//______________________________________________________________________________
void UMultiUnion::Extent(EAxisType aAxis, double& aMin, double& aMax) const
{
  // Determines the bounding box for the considered instance of "UMultipleUnion"
  UVector3 min, max;

  int numNodes = fSolids.size();
  for (int i = 0 ; i < numNodes ; ++i)
  {
    VUSolid& solid = *fSolids[i];
    UTransform3D transform = GetTransformation(i);
    solid.Extent(min, max);
   
    UUtils::TransformLimits(min, max, transform);
    
    if (i == 0)
    {
      switch (aAxis)
      {
        case eXaxis:
          aMin = min.x();
          aMax = max.x();
          break;
        case eYaxis:
          aMin = min.y();
          aMax = max.y();
          break;
        case eZaxis:
          aMin = min.z();
          aMax = max.z();
          break;
      }
    }
    else
    {
      // Deternine the min/max on the considered axis:
      switch (aAxis)
      {
        case eXaxis:
          if (min.x() < aMin)
            aMin = min.x();
          if (max.x() > aMax)
            aMax = max.x();
          break;
        case eYaxis:
          if (min.y() < aMin)
            aMin = min.y();
          if (max.y() > aMax)
            aMax = max.y();
          break;
        case eZaxis:
          if (min.z() < aMin)
            aMin = min.z();
          if (max.z() > aMax)
            aMax = max.z();
          break;
      }
    }
  }
}

//______________________________________________________________________________
void UMultiUnion::Extent(UVector3& aMin, UVector3& aMax) const
{
   Extent(eXaxis, aMin[0], aMax[0]);
   Extent(eYaxis, aMin[1], aMax[1]);
   Extent(eZaxis, aMin[2], aMax[2]);
}

//______________________________________________________________________________
bool UMultiUnion::Normal(const UVector3& aPoint, UVector3& aNormal) const
{
  // Computes the localNormal on a surface and returns it as a unit vector
  //   In case a point is further than toleranceNormal from a surface, set validNormal=false
  //   Must return a valid vector. (even if the point is not on the surface.)
  //
  //   On an edge or corner, provide an average localNormal of all facets within tolerance
  // NOTE: the tolerance value used in here is not yet the global surface
  //     tolerance - we will have to revise this value - TODO
  vector<int> candidates;
  UVector3 localPoint, normal, localNormal;
  double safety = UUtils::kInfinity;
  int node = -1;
  double normalTolerance = 1E-5;

  ///////////////////////////////////////////////////////////////////////////
  // Important comment: Cases for which the point is located on an edge or
  // on a vertice remain to be treated

  // determine weather we are in voxel area
  if (fVoxels.GetCandidatesVoxelArray(aPoint, candidates))
  {
    int limit = candidates.size();
    for (int i = 0 ; i < limit ; ++i)
    {
      int candidate = candidates[i];
      const UTransform3D& transform = fTransformObjs[candidate];
      // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
      localPoint = transform.LocalPoint(aPoint);
      VUSolid& solid = *fSolids[candidate];
      VUSolid::EnumInside location = solid.Inside(localPoint);

      if (location == eSurface)
      {
        // normal case when point is on surface, we pick first solid
        solid.Normal(localPoint, localNormal);
        normal = transform.GlobalVector(localNormal);
        aNormal = normal.Unit();
        return true;
      }
      else
      {
        // collect the smallest safety and remember solid node
        double s = (location == eInside) ? solid.SafetyFromInside(localPoint) : solid.SafetyFromOutside(localPoint);
        if (s < safety)
        {
          safety = s;
          node = candidate;
        }
      }
    }
    // on none of the solids, the point was not on the surface
    VUSolid& solid = *fSolids[node];
    const UTransform3D& transform = fTransformObjs[node];
    localPoint = transform.LocalPoint(aPoint);

    solid.Normal(localPoint, localNormal);
    normal = transform.GlobalVector(localNormal);
    aNormal = normal.Unit();
    if (safety > normalTolerance) return false;
    return true;
  }
  else
  {
    // for the case when point is certainly outside:

    // find a solid in union with the smallest safety
    node = SafetyFromOutsideNumberNode(aPoint, true, safety);
    VUSolid& solid = *fSolids[node];

    const UTransform3D& transform = fTransformObjs[node];
    localPoint = transform.LocalPoint(aPoint);

    // evaluate normal for point at this found solid
    solid.Normal(localPoint, localNormal);

    // transform multi-union coordinates
    normal = transform.GlobalVector(localNormal);

    aNormal = normal.Unit();
    if (safety > normalTolerance) return false;
    return true;
  }
}

//______________________________________________________________________________
double UMultiUnion::SafetyFromInside(const UVector3& point, bool aAccurate) const
{
  // Estimates isotropic distance to the surface of the solid. This must
  // be either accurate or an underestimate.
  //  Two modes: - default/fast mode, sacrificing accuracy for speed
  //             - "precise" mode,  requests accurate value if available.

  vector<int> candidates;
  UVector3 localPoint;
  double safetyMin = UUtils::kInfinity;

  // In general, the value return by SafetyFromInside will not be the exact
  // but only an undervalue (cf. overlaps)
  fVoxels.GetCandidatesVoxelArray(point, candidates);

  int limit = candidates.size();
  for (int i = 0; i < limit; ++i)
  {
    int candidate = candidates[i];
    // The coordinates of the point are modified so as to fit the intrinsic solid local frame:
    const UTransform3D& transform = fTransformObjs[candidate];
    localPoint = transform.LocalPoint(point);
    VUSolid& solid = *fSolids[candidate];
    if (solid.Inside(localPoint) == eInside)
    {
      double safety = solid.SafetyFromInside(localPoint, aAccurate);
      if (safetyMin > safety) safetyMin = safety;
    }
  }
  if (safetyMin == UUtils::kInfinity) safetyMin = 0; // we are not inside

  return safetyMin;
}

//______________________________________________________________________________
double UMultiUnion::SafetyFromOutside(const UVector3& point, bool aAccurate) const
{
  // Estimates the isotropic safety from a point outside the current solid to any
  // of its surfaces. The algorithm may be accurate or should provide a fast
  // underestimate.

  if (!aAccurate) return fVoxels.SafetyToBoundingBox(point);

  const std::vector<UVoxelBox>& boxes = fVoxels.GetBoxes();
  double safetyMin = UUtils::kInfinity;
  UVector3 localPoint;

  int numNodes = fSolids.size();
  for (int j = 0; j < numNodes; j++)
  {
    UVector3 dxyz;
    if (j > 0)
    {
      const UVector3& pos = boxes[j].pos;
      const UVector3& hlen = boxes[j].hlen;
      for (int i = 0; i <= 2; ++i)
        // distance to middle point - hlength => distance from point to border of x,y,z
        if ((dxyz[i] = std::abs(point[i] - pos[i]) - hlen[i]) > safetyMin)
          continue;

      double d2xyz = 0.;
      for (int i = 0; i <= 2; ++i)
        if (dxyz[i] > 0) d2xyz += dxyz[i] * dxyz[i];

      // minimal distance is at least this, but could be even higher. therefore, we can stop if previous was already lower, let us check if it does any chance to be better tha previous values...
      if (d2xyz >= safetyMin * safetyMin)
      {
#ifdef DEBUG
        const UTransform3D& transform = fTransformObjs[j];
        localPoint = transform.LocalPoint(point); // NOTE: ROOT does not make this transformation, although it does it at SafetyFromInside
        VUSolid& solid = *solids[j];
        double safety = solid.SafetyFromOutside(localPoint, true);
        if (safetyMin > safety)
          safety = safety;
#endif
        continue;
      }
    }
    const UTransform3D& transform = fTransformObjs[j];
    localPoint = transform.LocalPoint(point); // NOTE: ROOT does not make this transformation, although it does it at SafetyFromInside
    VUSolid& solid = *fSolids[j];

    double safety = solid.SafetyFromOutside(localPoint, aAccurate); // careful, with aAcurate it can return underestimate, than the condition d2xyz >= safetyMin*safetyMin does not return same result as Geant4 or Root boolean union, it actually return better values, more close to surface. so although the values of this method would be correct, although different from Geant4
    if (safety <= 0) return safety; // it was detected, that the point is not located outside
    if (safetyMin > safety) safetyMin = safety;
  }
  return safetyMin;
}

//______________________________________________________________________________
double UMultiUnion::SurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    ;
  }
  else
  {
    fSurfaceArea = EstimateSurfaceArea(1000000, 0.001);
  }
  return fSurfaceArea;

}

//______________________________________________________________________________
void UMultiUnion::Voxelize()
{
  fVoxels.Voxelize(fSolids, fTransformObjs);
}

//______________________________________________________________________________
int UMultiUnion::SafetyFromOutsideNumberNode(const UVector3& aPoint, bool /*aAccurate*/, double& safetyMin) const
{
  // Method returning the closest node from a point located outside a UMultiUnion.
  // This is used to compute the normal in the case no candidate has been found.

  const std::vector<UVoxelBox>& boxes = fVoxels.GetBoxes();
  safetyMin = UUtils::kInfinity;
  int safetyNode = -1;
  UVector3 localPoint;

  int numNodes = fSolids.size();
  for (int i = 0; i < numNodes; ++i)
  {
    double d2xyz = 0.;
    double dxyz0 = std::abs(aPoint.x() - boxes[i].pos.x()) - boxes[i].hlen.x();
    if (dxyz0 > safetyMin) continue;
    double dxyz1 = std::abs(aPoint.y() - boxes[i].pos.y()) - boxes[i].hlen.y();
    if (dxyz1 > safetyMin) continue;
    double dxyz2 = std::abs(aPoint.z() - boxes[i].pos.z()) - boxes[i].hlen.z();
    if (dxyz2 > safetyMin) continue;

    if (dxyz0 > 0) d2xyz += dxyz0 * dxyz0;
    if (dxyz1 > 0) d2xyz += dxyz1 * dxyz1;
    if (dxyz2 > 0) d2xyz += dxyz2 * dxyz2;
    if (d2xyz >= safetyMin * safetyMin) continue;

    VUSolid& solid = *fSolids[i];
    const UTransform3D& transform = fTransformObjs[i];
    localPoint = transform.LocalPoint(aPoint);
    double safety = solid.SafetyFromOutside(localPoint, true);
    if (safetyMin > safety)
    {
      safetyMin = safety;
      safetyNode = i;
    }
  }
  return safetyNode;
}

// Stream object contents to an output stream
//______________________________________________________________________________
std::ostream& UMultiUnion::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "                *** Dump for solid - " << GetName() << " ***\n"
     << "                ===================================================\n"
     << " Solid type: UMultiUnion\n"
     << " Parameters: \n";
     int numNodes = fSolids.size();
     for (int i = 0 ; i < numNodes ; ++i)
     {
      VUSolid& solid = *fSolids[i];
      solid.StreamInfo(os);
      const UTransform3D& transform = fTransformObjs[i];
      os<<" Translation is " <<transform.fTr<<" \n";
      os<<" Rotation is :"<<" \n";
      for(int j = 0; j < 3; j++ )
	os<<" "<<transform.fRot[j*3]<<" "<<transform.fRot[1+j*3]<<" "<<transform.fRot[2+j*3]<<"\n";
     }

  os << "             \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//______________________________________________________________________________
UVector3 UMultiUnion::GetPointOnSurface() const
{
  UVector3 point;

  int size = fSolids.size();

  do
  {
    int random = (int) UUtils::Random(0, size);
    /*
    for (int i = 0; i < size; i++)
    {
      VUSolid &solid = *fSolids[i];
    }
    */
    VUSolid& solid = *fSolids[random];
    point = solid.GetPointOnSurface();
    const UTransform3D& transform = fTransformObjs[random];
    point = transform.GlobalPoint(point);
  }
  while (Inside(point) != eSurface);

  return point;
}
