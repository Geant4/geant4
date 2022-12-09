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
// Implementation of G4BoundingEnvelope
//
// 2016.05.25, E.Tcherniaev - initial version
// --------------------------------------------------------------------

#include <cmath>

#include "globals.hh"
#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"
#include "G4Normal3D.hh"

const G4double kCarTolerance =
  G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

///////////////////////////////////////////////////////////////////////
//
// Constructor from an axis aligned bounding box
//
G4BoundingEnvelope::G4BoundingEnvelope(const G4ThreeVector& pMin,
                                       const G4ThreeVector& pMax)
  : fMin(pMin), fMax(pMax)
{
  // Check correctness of bounding box
  //
  CheckBoundingBox();
}

///////////////////////////////////////////////////////////////////////
//
// Constructor from a sequence of polygons
//
G4BoundingEnvelope::
G4BoundingEnvelope(const std::vector<const G4ThreeVectorList*>& polygons)
  : fPolygons(&polygons)
{
  // Check correctness of polygons
  //
  CheckBoundingPolygons();

  // Set bounding box
  //
  G4double xmin =  kInfinity, ymin =  kInfinity, zmin =  kInfinity;
  G4double xmax = -kInfinity, ymax = -kInfinity, zmax = -kInfinity;
  for (auto ibase = fPolygons->cbegin(); ibase != fPolygons->cend(); ++ibase)
  {
    for (auto ipoint = (*ibase)->cbegin(); ipoint != (*ibase)->cend(); ++ipoint)
    {
      G4double x = ipoint->x();
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      G4double y = ipoint->y();
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
      G4double z = ipoint->z();
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
    }
  }
  fMin.set(xmin,ymin,zmin);
  fMax.set(xmax,ymax,zmax);

  // Check correctness of bounding box
  //
  CheckBoundingBox();
}

///////////////////////////////////////////////////////////////////////
//
// Constructor from a bounding box and a sequence of polygons
//
G4BoundingEnvelope::
G4BoundingEnvelope( const G4ThreeVector& pMin,
                    const G4ThreeVector& pMax,
                    const std::vector<const G4ThreeVectorList*>& polygons)
  : fMin(pMin), fMax(pMax), fPolygons(&polygons)
{
  // Check correctness of bounding box and polygons
  //
  CheckBoundingBox();
  CheckBoundingPolygons();
}

///////////////////////////////////////////////////////////////////////
//
// Check correctness of the axis aligned bounding box
//
void G4BoundingEnvelope::CheckBoundingBox()
{
  if (fMin.x() >= fMax.x() || fMin.y() >= fMax.y() || fMin.z() >= fMax.z())
  {
    std::ostringstream message;
    message << "Badly defined bounding box (min >= max)!"
            << "\npMin = " << fMin
            << "\npMax = " << fMax;
    G4Exception("G4BoundingEnvelope::CheckBoundingBox()",
                "GeomMgt0001", JustWarning, message);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Check correctness of the sequence of bounding polygons.
// Firsf and last polygons may consist of a single vertex
//
void G4BoundingEnvelope::CheckBoundingPolygons()
{
  std::size_t nbases = fPolygons->size();
  if (nbases < 2)
  {
    std::ostringstream message;
    message << "Wrong number of polygons in the sequence: " << nbases
            << "\nShould be at least two!";
    G4Exception("G4BoundingEnvelope::CheckBoundingPolygons()",
                "GeomMgt0001", FatalException, message);
    return;
  }

  std::size_t nsize  = std::max((*fPolygons)[0]->size(),(*fPolygons)[1]->size());
  if (nsize < 3)
  {
    std::ostringstream message;
    message << "Badly constructed polygons!"
            << "\nNumber of polygons: " << nbases
            << "\nPolygon #0 size: " << (*fPolygons)[0]->size()
            << "\nPolygon #1 size: " << (*fPolygons)[1]->size()
            << "\n...";
    G4Exception("G4BoundingEnvelope::CheckBoundingPolygons()",
                "GeomMgt0001", FatalException, message);
    return;
  }

  for (std::size_t k=0; k<nbases; ++k)
  {
    std::size_t np = (*fPolygons)[k]->size();
    if (np == nsize)            continue;
    if (np == 1 && k==0)        continue;
    if (np == 1 && k==nbases-1) continue;
    std::ostringstream message;
    message << "Badly constructed polygons!"
            << "\nNumber of polygons: " << nbases
            << "\nPolygon #" << k << " size: " << np
            << "\nexpected size: " << nsize;
    G4Exception("G4BoundingEnvelope::SetBoundingPolygons()",
                "GeomMgt0001", FatalException, message);
    return;
  }
}

///////////////////////////////////////////////////////////////////////
//
// Quick comparison: bounding box vs voxel, it return true if further
// calculations are not needed
//
G4bool
G4BoundingEnvelope::
BoundingBoxVsVoxelLimits(const EAxis pAxis,
                         const G4VoxelLimits& pVoxelLimits,
                         const G4Transform3D& pTransform3D,
                         G4double& pMin, G4double& pMax) const
{
  pMin =  kInfinity;
  pMax = -kInfinity;
  G4double xminlim = pVoxelLimits.GetMinXExtent();
  G4double xmaxlim = pVoxelLimits.GetMaxXExtent();
  G4double yminlim = pVoxelLimits.GetMinYExtent();
  G4double ymaxlim = pVoxelLimits.GetMaxYExtent();
  G4double zminlim = pVoxelLimits.GetMinZExtent();
  G4double zmaxlim = pVoxelLimits.GetMaxZExtent();

  // Special case of pure translation
  //
  if (pTransform3D.xx()==1 && pTransform3D.yy()==1 && pTransform3D.zz()==1)
  {
    G4double xmin = fMin.x() + pTransform3D.dx();
    G4double xmax = fMax.x() + pTransform3D.dx();
    G4double ymin = fMin.y() + pTransform3D.dy();
    G4double ymax = fMax.y() + pTransform3D.dy();
    G4double zmin = fMin.z() + pTransform3D.dz();
    G4double zmax = fMax.z() + pTransform3D.dz();

    if (xmin-kCarTolerance > xmaxlim) return true;
    if (xmax+kCarTolerance < xminlim) return true;
    if (ymin-kCarTolerance > ymaxlim) return true;
    if (ymax+kCarTolerance < yminlim) return true;
    if (zmin-kCarTolerance > zmaxlim) return true;
    if (zmax+kCarTolerance < zminlim) return true;

    if (xmin >= xminlim && xmax <= xmaxlim &&
        ymin >= yminlim && ymax <= ymaxlim &&
        zmin >= zminlim && zmax <= zmaxlim)
    {
      if (pAxis == kXAxis)
      {
        pMin = (xmin-kCarTolerance < xminlim) ? xminlim : xmin;
        pMax = (xmax+kCarTolerance > xmaxlim) ? xmaxlim : xmax;
      }
      else if (pAxis == kYAxis)
      {
        pMin = (ymin-kCarTolerance < yminlim) ? yminlim : ymin;
        pMax = (ymax+kCarTolerance > ymaxlim) ? ymaxlim : ymax;
      }
      else if (pAxis == kZAxis)
      {
        pMin = (zmin-kCarTolerance < zminlim) ? zminlim : zmin;
        pMax = (zmax+kCarTolerance > zmaxlim) ? zmaxlim : zmax;
      }
      pMin -= kCarTolerance;
      pMax += kCarTolerance;
      return true;
    }
  }

  // Find max scale factor of the transformation, set delta
  // equal to kCarTolerance multiplied by the scale factor
  //
  G4double scale = FindScaleFactor(pTransform3D);
  G4double delta = kCarTolerance*scale;

  // Set the sphere surrounding the bounding box
  //
  G4Point3D center = pTransform3D*G4Point3D(0.5*(fMin+fMax));
  G4double  radius = 0.5*scale*(fMax-fMin).mag() + delta;

  // Check if the sphere surrounding the bounding box is outside
  // the voxel limits
  //
  if (center.x()-radius > xmaxlim) return true;
  if (center.y()-radius > ymaxlim) return true;
  if (center.z()-radius > zmaxlim) return true;
  if (center.x()+radius < xminlim) return true;
  if (center.y()+radius < yminlim) return true;
  if (center.z()+radius < zminlim) return true;
  return false;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate extent of the specified bounding envelope
//
G4bool
G4BoundingEnvelope::CalculateExtent(const EAxis pAxis,
                                    const G4VoxelLimits& pVoxelLimits,
                                    const G4Transform3D& pTransform3D,
                                    G4double& pMin, G4double& pMax) const
{
  pMin =  kInfinity;
  pMax = -kInfinity;
  G4double xminlim = pVoxelLimits.GetMinXExtent();
  G4double xmaxlim = pVoxelLimits.GetMaxXExtent();
  G4double yminlim = pVoxelLimits.GetMinYExtent();
  G4double ymaxlim = pVoxelLimits.GetMaxYExtent();
  G4double zminlim = pVoxelLimits.GetMinZExtent();
  G4double zmaxlim = pVoxelLimits.GetMaxZExtent();

  // Special case of pure translation
  //
  if (pTransform3D.xx()==1 && pTransform3D.yy()==1 && pTransform3D.zz()==1)
  {
    G4double xmin = fMin.x() + pTransform3D.dx();
    G4double xmax = fMax.x() + pTransform3D.dx();
    G4double ymin = fMin.y() + pTransform3D.dy();
    G4double ymax = fMax.y() + pTransform3D.dy();
    G4double zmin = fMin.z() + pTransform3D.dz();
    G4double zmax = fMax.z() + pTransform3D.dz();

    if (xmin-kCarTolerance > xmaxlim) return false;
    if (xmax+kCarTolerance < xminlim) return false;
    if (ymin-kCarTolerance > ymaxlim) return false;
    if (ymax+kCarTolerance < yminlim) return false;
    if (zmin-kCarTolerance > zmaxlim) return false;
    if (zmax+kCarTolerance < zminlim) return false;

    if (fPolygons == nullptr)
    {
      if (pAxis == kXAxis)
      {
        pMin = (xmin-kCarTolerance < xminlim) ? xminlim : xmin;
        pMax = (xmax+kCarTolerance > xmaxlim) ? xmaxlim : xmax;
      }
      else if (pAxis == kYAxis)
      {
        pMin = (ymin-kCarTolerance < yminlim) ? yminlim : ymin;
        pMax = (ymax+kCarTolerance > ymaxlim) ? ymaxlim : ymax;
      }
      else if (pAxis == kZAxis)
      {
        pMin = (zmin-kCarTolerance < zminlim) ? zminlim : zmin;
        pMax = (zmax+kCarTolerance > zmaxlim) ? zmaxlim : zmax;
      }
      pMin -= kCarTolerance;
      pMax += kCarTolerance;
      return true;
    }
  }

  // Find max scale factor of the transformation, set delta
  // equal to kCarTolerance multiplied by the scale factor
  //
  G4double scale = FindScaleFactor(pTransform3D);
  G4double delta = kCarTolerance*scale;

  // Set the sphere surrounding the bounding box
  //
  G4Point3D center = pTransform3D*G4Point3D(0.5*(fMin+fMax));
  G4double  radius = 0.5*scale*(fMax-fMin).mag() + delta;

  // Check if the sphere surrounding the bounding box is within
  // the voxel limits, if so then transform only one coordinate
  //
  if (center.x()-radius >= xminlim && center.x()+radius <= xmaxlim &&
      center.y()-radius >= yminlim && center.y()+radius <= ymaxlim &&
      center.z()-radius >= zminlim && center.z()+radius <= zmaxlim )
  {
    G4double cx, cy, cz, cd;
    if (pAxis == kXAxis)
    {
      cx = pTransform3D.xx();
      cy = pTransform3D.xy();
      cz = pTransform3D.xz();
      cd = pTransform3D.dx();
    }
    else if (pAxis == kYAxis)
    {
      cx = pTransform3D.yx();
      cy = pTransform3D.yy();
      cz = pTransform3D.yz();
      cd = pTransform3D.dy();
    }
    else if (pAxis == kZAxis)
    {
      cx = pTransform3D.zx();
      cy = pTransform3D.zy();
      cz = pTransform3D.zz();
      cd = pTransform3D.dz();
    }
    else
    {
      cx = cy = cz = cd = kInfinity;
    }
    G4double emin = kInfinity, emax = -kInfinity;
    if (fPolygons == nullptr)
    {
      G4double coor;
      coor = cx*fMin.x() + cy*fMin.y() + cz*fMin.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMax.x() + cy*fMin.y() + cz*fMin.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMax.x() + cy*fMax.y() + cz*fMin.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMin.x() + cy*fMax.y() + cz*fMin.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMin.x() + cy*fMin.y() + cz*fMax.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMax.x() + cy*fMin.y() + cz*fMax.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMax.x() + cy*fMax.y() + cz*fMax.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
      coor = cx*fMin.x() + cy*fMax.y() + cz*fMax.z() + cd;
      if (coor < emin) emin = coor;
      if (coor > emax) emax = coor;
    }
    else
    {
      for (auto ibase=fPolygons->cbegin(); ibase!=fPolygons->cend(); ++ibase)
      {
        for (auto ipoint=(*ibase)->cbegin(); ipoint!=(*ibase)->cend(); ++ipoint)
        {
          G4double coor = ipoint->x()*cx + ipoint->y()*cy + ipoint->z()*cz + cd;
          if (coor < emin) emin = coor;
          if (coor > emax) emax = coor;
        }
      }
    }
    pMin = emin - delta;
    pMax = emax + delta;
    return true;
  }

  // Check if the sphere surrounding the bounding box is outside
  // the voxel limits
  //
  if (center.x()-radius > xmaxlim) return false;
  if (center.y()-radius > ymaxlim) return false;
  if (center.z()-radius > zmaxlim) return false;
  if (center.x()+radius < xminlim) return false;
  if (center.y()+radius < yminlim) return false;
  if (center.z()+radius < zminlim) return false;

  // Transform polygons
  //
  std::vector<G4Point3D> vertices;
  std::vector<std::pair<G4int, G4int>> bases;
  TransformVertices(pTransform3D, vertices, bases);
  std::size_t nbases = bases.size();

  // Create adjusted G4VoxelLimits box. New limits are extended by
  // delta, kCarTolerance multiplied by max scale factor of
  // the transformation
  //
  EAxis axis[] = { kXAxis, kYAxis, kZAxis };
  G4VoxelLimits limits; // default is unlimited
  for (auto i=0; i<3; ++i)
  {
    if (pVoxelLimits.IsLimited(axis[i]))
    {
      G4double emin = pVoxelLimits.GetMinExtent(axis[i]) - delta;
      G4double emax = pVoxelLimits.GetMaxExtent(axis[i]) + delta;
      limits.AddLimit(axis[i], emin, emax);
    }
  }

  // Main loop along the set of prisms
  //
  G4Polygon3D baseA, baseB;
  G4Segment3D extent;
  extent.first  = G4Point3D( kInfinity, kInfinity, kInfinity);
  extent.second = G4Point3D(-kInfinity,-kInfinity,-kInfinity);
  for (std::size_t k=0; k<nbases-1; ++k)
  {
    baseA.resize(bases[k].second);
    for (G4int i = 0; i < bases[k].second; ++i)
      baseA[i] = vertices[bases[k].first + i];

    baseB.resize(bases[k+1].second);
    for (G4int i = 0; i < bases[k+1].second; ++i)
      baseB[i] = vertices[bases[k+1].first + i];

    // Find bounding box of current prism
    G4Segment3D  prismAABB;
    GetPrismAABB(baseA, baseB, prismAABB);

    // Check if prismAABB is completely within the voxel limits
    if (prismAABB.first.x() >= limits.GetMinXExtent() &&
        prismAABB.first.y() >= limits.GetMinYExtent() &&
        prismAABB.first.z() >= limits.GetMinZExtent() &&
        prismAABB.second.x()<= limits.GetMaxXExtent() &&
        prismAABB.second.y()<= limits.GetMaxYExtent() &&
        prismAABB.second.z()<= limits.GetMaxZExtent())
    {
      if (extent.first.x()  > prismAABB.first.x())
        extent.first.setX( prismAABB.first.x() );
      if (extent.first.y()  > prismAABB.first.y())
        extent.first.setY( prismAABB.first.y() );
      if (extent.first.z()  > prismAABB.first.z())
        extent.first.setZ( prismAABB.first.z() );
      if (extent.second.x() < prismAABB.second.x())
        extent.second.setX(prismAABB.second.x());
      if (extent.second.y() < prismAABB.second.y())
        extent.second.setY(prismAABB.second.y());
      if (extent.second.z() < prismAABB.second.z())
        extent.second.setZ(prismAABB.second.z());
      continue;
    }

    // Check if prismAABB is outside the voxel limits
    if (prismAABB.first.x()  > limits.GetMaxXExtent()) continue;
    if (prismAABB.first.y()  > limits.GetMaxYExtent()) continue;
    if (prismAABB.first.z()  > limits.GetMaxZExtent()) continue;
    if (prismAABB.second.x() < limits.GetMinXExtent()) continue;
    if (prismAABB.second.y() < limits.GetMinYExtent()) continue;
    if (prismAABB.second.z() < limits.GetMinZExtent()) continue;

    // Clip edges of the prism by adjusted G4VoxelLimits box
    std::vector<G4Segment3D> vecEdges;
    CreateListOfEdges(baseA, baseB, vecEdges);
    if (ClipEdgesByVoxel(vecEdges, limits, extent)) continue;

    // Some edges of the prism are completely outside of the voxel
    // limits, clip selected edges (see bits) of adjusted G4VoxelLimits
    // by the prism
    G4int bits = 0x000;
    if (limits.GetMinXExtent() < prismAABB.first.x())
      bits |= 0x988; // 1001 1000 1000
    if (limits.GetMaxXExtent() > prismAABB.second.x())
      bits |= 0x622; // 0110 0010 0010

    if (limits.GetMinYExtent() < prismAABB.first.y())
      bits |= 0x311; // 0011 0001 0001
    if (limits.GetMaxYExtent() > prismAABB.second.y())
      bits |= 0xC44; // 1100 0100 0100

    if (limits.GetMinZExtent() < prismAABB.first.z())
      bits |= 0x00F; // 0000 0000 1111
    if (limits.GetMaxZExtent() > prismAABB.second.z())
      bits |= 0x0F0; // 0000 1111 0000
    if (bits == 0xFFF) continue;

    std::vector<G4Plane3D> vecPlanes;
    CreateListOfPlanes(baseA, baseB, vecPlanes);
    ClipVoxelByPlanes(bits, limits, vecPlanes, prismAABB, extent);
  } // End of the main loop

  // Final adjustment of the extent
  //
  G4double emin = 0, emax = 0;
  if (pAxis == kXAxis) { emin = extent.first.x(); emax = extent.second.x(); }
  if (pAxis == kYAxis) { emin = extent.first.y(); emax = extent.second.y(); }
  if (pAxis == kZAxis) { emin = extent.first.z(); emax = extent.second.z(); }

  if (emin > emax) return false;
  emin -= delta;
  emax += delta;
  G4double minlim = pVoxelLimits.GetMinExtent(pAxis);
  G4double maxlim = pVoxelLimits.GetMaxExtent(pAxis);
  pMin = (emin < minlim) ? minlim-kCarTolerance : emin;
  pMax = (emax > maxlim) ? maxlim+kCarTolerance : emax;
  return true;
}

///////////////////////////////////////////////////////////////////////
//
// Find max scale factor of the transformation
//
G4double
G4BoundingEnvelope::FindScaleFactor(const G4Transform3D& pTransform3D) const
{
  if (pTransform3D.xx() == 1. &&
      pTransform3D.yy() == 1. &&
      pTransform3D.zz() == 1.) return 1.;

  G4double xx = pTransform3D.xx();
  G4double yx = pTransform3D.yx();
  G4double zx = pTransform3D.zx();
  G4double sxsx = xx*xx + yx*yx + zx*zx;
  G4double xy = pTransform3D.xy();
  G4double yy = pTransform3D.yy();
  G4double zy = pTransform3D.zy();
  G4double sysy = xy*xy + yy*yy + zy*zy;
  G4double xz = pTransform3D.xz();
  G4double yz = pTransform3D.yz();
  G4double zz = pTransform3D.zz();
  G4double szsz = xz*xz + yz*yz + zz*zz;
  G4double ss = std::max(std::max(sxsx,sysy),szsz);
  return (ss <= 1.) ? 1. : std::sqrt(ss);
}

///////////////////////////////////////////////////////////////////////
//
// Transform polygonal bases
//
void
G4BoundingEnvelope::
TransformVertices(const G4Transform3D& pTransform3D,
                  std::vector<G4Point3D>& pVertices,
                  std::vector<std::pair<G4int, G4int>>& pBases) const
{
  G4ThreeVectorList baseA(4), baseB(4);
  std::vector<const G4ThreeVectorList*> aabb(2);
  aabb[0] = &baseA;
  aabb[1] = &baseB;
  if (fPolygons == nullptr)
  {
    baseA[0].set(fMin.x(),fMin.y(),fMin.z());
    baseA[1].set(fMax.x(),fMin.y(),fMin.z());
    baseA[2].set(fMax.x(),fMax.y(),fMin.z());
    baseA[3].set(fMin.x(),fMax.y(),fMin.z());
    baseB[0].set(fMin.x(),fMin.y(),fMax.z());
    baseB[1].set(fMax.x(),fMin.y(),fMax.z());
    baseB[2].set(fMax.x(),fMax.y(),fMax.z());
    baseB[3].set(fMin.x(),fMax.y(),fMax.z());
  }
  auto ia    = (fPolygons == nullptr) ? aabb.cbegin() : fPolygons->cbegin();
  auto iaend = (fPolygons == nullptr) ? aabb.cend()   : fPolygons->cend();

  // Fill vector of bases
  //
  G4int index = 0;
  for (auto i = ia; i != iaend; ++i)
  {
    G4int nv = (G4int)(*i)->size();
    pBases.push_back(std::make_pair(index, nv));
    index += nv;
  }

  // Fill vector of transformed vertices
  //
  if (pTransform3D.xx() == 1. &&
      pTransform3D.yy() == 1. &&
      pTransform3D.zz() == 1.)
  {
    G4ThreeVector offset = pTransform3D.getTranslation();
    for (auto i = ia; i != iaend; ++i)
      for (auto k = (*i)->cbegin(); k != (*i)->cend(); ++k)
        pVertices.push_back(G4Point3D((*k) + offset));
  }
  else
  {
    for (auto i = ia; i != iaend; ++i)
      for (auto k = (*i)->cbegin(); k != (*i)->cend(); ++k)
        pVertices.push_back(pTransform3D*G4Point3D(*k));
  }
}

///////////////////////////////////////////////////////////////////////
//
// Find bounding box of a prism
//
void
G4BoundingEnvelope::GetPrismAABB(const G4Polygon3D& pBaseA,
                                 const G4Polygon3D& pBaseB,
                                       G4Segment3D& pAABB) const
{
  G4double xmin =  kInfinity, ymin =  kInfinity, zmin =  kInfinity;
  G4double xmax = -kInfinity, ymax = -kInfinity, zmax = -kInfinity;

  // First base
  //
  for (auto it1 = pBaseA.cbegin(); it1 != pBaseA.cend(); ++it1)
  {
    G4double x = it1->x();
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    G4double y = it1->y();
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    G4double z = it1->z();
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }

  // Second base
  //
  for (auto it2 = pBaseB.cbegin(); it2 != pBaseB.cend(); ++it2)
  {
    G4double x = it2->x();
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    G4double y = it2->y();
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    G4double z = it2->z();
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
  }

  // Set bounding box
  //
  pAABB.first  = G4Point3D(xmin,ymin,zmin);
  pAABB.second = G4Point3D(xmax,ymax,zmax);
}

///////////////////////////////////////////////////////////////////////
//
// Create list of edges of a prism
//
void
G4BoundingEnvelope::CreateListOfEdges(const G4Polygon3D& baseA,
                                      const G4Polygon3D& baseB,
                                      std::vector<G4Segment3D>& pEdges) const
{
  std::size_t na = baseA.size();
  std::size_t nb = baseB.size();
  pEdges.clear();
  if (na == nb)
  {
    pEdges.resize(3*na);
    std::size_t k = na - 1;
    for (std::size_t i=0; i<na; ++i)
    {
      pEdges.push_back(G4Segment3D(baseA[i],baseB[i]));
      pEdges.push_back(G4Segment3D(baseA[i],baseA[k]));
      pEdges.push_back(G4Segment3D(baseB[i],baseB[k]));
      k = i;
    }
  }
  else if (nb == 1)
  {
    pEdges.resize(2*na);
    std::size_t k = na - 1;
    for (std::size_t i=0; i<na; ++i)
    {
      pEdges.push_back(G4Segment3D(baseA[i],baseA[k]));
      pEdges.push_back(G4Segment3D(baseA[i],baseB[0]));
      k = i;
    }
  }
  else if (na == 1)
  {
    pEdges.resize(2*nb);
    std::size_t k = nb - 1;
    for (std::size_t i=0; i<nb; ++i)
    {
      pEdges.push_back(G4Segment3D(baseB[i],baseB[k]));
      pEdges.push_back(G4Segment3D(baseB[i],baseA[0]));
      k = i;
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Create list of planes bounding a prism
//
void
G4BoundingEnvelope::CreateListOfPlanes(const G4Polygon3D& baseA,
                                       const G4Polygon3D& baseB,
                                       std::vector<G4Plane3D>& pPlanes) const
{
  // Find centers of the bases and internal point of the prism
  //
  std::size_t na = baseA.size();
  std::size_t nb = baseB.size();
  G4Point3D pa(0.,0.,0.), pb(0.,0.,0.), p0;
  G4Normal3D norm;
  for (std::size_t i=0; i<na; ++i) pa += baseA[i];
  for (std::size_t i=0; i<nb; ++i) pb += baseB[i];
  pa /= na; pb /= nb; p0 = (pa+pb)/2.;

  // Create list of planes
  //
  pPlanes.clear();
  if (na == nb)  // bases with equal number of vertices
  {
    std::size_t k = na - 1;
    for (std::size_t i=0; i<na; ++i)
    {
      norm = (baseB[k]-baseA[i]).cross(baseA[k]-baseB[i]);
      if (norm.mag2() > kCarTolerance)
      {
        pPlanes.push_back(G4Plane3D(norm,baseA[i]));
      }
      k = i;
    }
    norm = (baseA[2]-baseA[0]).cross(baseA[1]-pa);
    if (norm.mag2() > kCarTolerance)
    {
      pPlanes.push_back(G4Plane3D(norm,pa));
    }
    norm = (baseB[2]-baseB[0]).cross(baseB[1]-pb);
    if (norm.mag2() > kCarTolerance)
    {
      pPlanes.push_back(G4Plane3D(norm,pb));
    }
  }
  else if (nb == 1) // baseB has one vertex
  {
    std::size_t k = na - 1;
    for (std::size_t i=0; i<na; ++i)
    {
      norm = (baseA[i]-baseB[0]).cross(baseA[k]-baseB[0]);
      if (norm.mag2() > kCarTolerance)
      {
        pPlanes.push_back(G4Plane3D(norm,baseB[0]));
      }
      k = i;
    }
    norm = (baseA[2]-baseA[0]).cross(baseA[1]-pa);
    if (norm.mag2() > kCarTolerance)
    {
      pPlanes.push_back(G4Plane3D(norm,pa));
    }
  }
  else if (na == 1) // baseA has one vertex
  {
    std::size_t k = nb - 1;
    for (std::size_t i=0; i<nb; ++i)
    {
      norm = (baseB[i]-baseA[0]).cross(baseB[k]-baseA[0]);
      if (norm.mag2() > kCarTolerance)
      {
        pPlanes.push_back(G4Plane3D(norm,baseA[0]));
      }
      k = i;
    }
    norm = (baseB[2]-baseB[0]).cross(baseB[1]-pb);
    if (norm.mag2() > kCarTolerance)
    {
      pPlanes.push_back(G4Plane3D(norm,pb));
    }
  }

  // Ensure that normals of the planes point to outside
  //
  std::size_t nplanes = pPlanes.size();
  for (std::size_t i=0; i<nplanes; ++i)
  {
    pPlanes[i].normalize();
    if (pPlanes[i].distance(p0) > 0)
    {
      pPlanes[i] = G4Plane3D(-pPlanes[i].a(),-pPlanes[i].b(),
                             -pPlanes[i].c(),-pPlanes[i].d());
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Clip edges of a prism by G4VoxelLimits box. Return true if all edges
// are inside or intersect the voxel, in this case further calculations
// are not needed
//
G4bool
G4BoundingEnvelope::ClipEdgesByVoxel(const std::vector<G4Segment3D>& pEdges,
                                     const G4VoxelLimits& pBox,
                                           G4Segment3D& pExtent) const
{
  G4bool    done = true;
  G4Point3D emin = pExtent.first;
  G4Point3D emax = pExtent.second;

  std::size_t nedges = pEdges.size();
  for (std::size_t k=0; k<nedges; ++k)
  {
    G4Point3D p1 = pEdges[k].first;
    G4Point3D p2 = pEdges[k].second;
    if (std::abs(p1.x()-p2.x())+
        std::abs(p1.y()-p2.y())+
        std::abs(p1.z()-p2.z()) < kCarTolerance) continue;
    G4double  d1, d2;
    // Clip current edge by X min
    d1 = pBox.GetMinXExtent() - p1.x();
    d2 = pBox.GetMinXExtent() - p2.x();
    if (d1 > 0.0)
    {
      if (d2 > 0.0) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);                   // move p1
    }
    else
    {
      if (d2 > 0.0) { p2 = (p1*d2-p2*d1)/(d2-d1); } // move p2
    }

    // Clip current edge by X max
    d1 = p1.x() - pBox.GetMaxXExtent();
    d2 = p2.x() - pBox.GetMaxXExtent();
    if (d1 > 0.)
    {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    }
    else
    {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Y min
    d1 = pBox.GetMinYExtent() - p1.y();
    d2 = pBox.GetMinYExtent() - p2.y();
    if (d1 > 0.)
    {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    }
    else
    {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Y max
    d1 = p1.y() - pBox.GetMaxYExtent();
    d2 = p2.y() - pBox.GetMaxYExtent();
    if (d1 > 0.)
    {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    }
    else
    {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Z min
    d1 = pBox.GetMinZExtent() - p1.z();
    d2 = pBox.GetMinZExtent() - p2.z();
    if (d1 > 0.)
    {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    }
    else
    {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Z max
    d1 = p1.z() - pBox.GetMaxZExtent();
    d2 = p2.z() - pBox.GetMaxZExtent();
    if (d1 > 0.)
    {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    }
    else
    {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Adjust current extent
    emin.setX(std::min(std::min(p1.x(),p2.x()),emin.x()));
    emin.setY(std::min(std::min(p1.y(),p2.y()),emin.y()));
    emin.setZ(std::min(std::min(p1.z(),p2.z()),emin.z()));

    emax.setX(std::max(std::max(p1.x(),p2.x()),emax.x()));
    emax.setY(std::max(std::max(p1.y(),p2.y()),emax.y()));
    emax.setZ(std::max(std::max(p1.z(),p2.z()),emax.z()));
  }

  // Return true if all edges (at least partially) are inside
  // the voxel limits, otherwise return false
  pExtent.first  = emin;
  pExtent.second = emax;

  return done;
}

///////////////////////////////////////////////////////////////////////
//
// Clip G4VoxelLimits by set of planes bounding a prism
//
void
G4BoundingEnvelope::ClipVoxelByPlanes(G4int pBits,
                                      const G4VoxelLimits& pBox,
                                      const std::vector<G4Plane3D>& pPlanes,
                                      const G4Segment3D& pAABB,
                                            G4Segment3D& pExtent) const
{
  G4Point3D emin = pExtent.first;
  G4Point3D emax = pExtent.second;

  // Create edges of the voxel limits box reducing them where
  // appropriate to avoid calculations with big numbers (kInfinity)
  //
  G4double xmin = std::max(pBox.GetMinXExtent(),pAABB.first.x() -1.);
  G4double xmax = std::min(pBox.GetMaxXExtent(),pAABB.second.x()+1.);

  G4double ymin = std::max(pBox.GetMinYExtent(),pAABB.first.y() -1.);
  G4double ymax = std::min(pBox.GetMaxYExtent(),pAABB.second.y()+1.);

  G4double zmin = std::max(pBox.GetMinZExtent(),pAABB.first.z() -1.);
  G4double zmax = std::min(pBox.GetMaxZExtent(),pAABB.second.z()+1.);

  std::vector<G4Segment3D> edges(12);
  G4int i = 0, bits = pBits;
  if (!(bits & 0x001))
  {
    edges[i  ].first.set( xmin,ymin,zmin);
    edges[i++].second.set(xmax,ymin,zmin);
  }
  if (!(bits & 0x002))
  {
    edges[i  ].first.set( xmax,ymin,zmin);
    edges[i++].second.set(xmax,ymax,zmin);
  }
  if (!(bits & 0x004))
  {
    edges[i  ].first.set( xmax,ymax,zmin);
    edges[i++].second.set(xmin,ymax,zmin);
  }
  if (!(bits & 0x008))
  {
    edges[i  ].first.set( xmin,ymax,zmin);
    edges[i++].second.set(xmin,ymin,zmin);
  }

  if (!(bits & 0x010))
  {
    edges[i  ].first.set( xmin,ymin,zmax);
    edges[i++].second.set(xmax,ymin,zmax);
  }
  if (!(bits & 0x020))
  {
    edges[i  ].first.set( xmax,ymin,zmax);
    edges[i++].second.set(xmax,ymax,zmax);
  }
  if (!(bits & 0x040))
  {
    edges[i  ].first.set( xmax,ymax,zmax);
    edges[i++].second.set(xmin,ymax,zmax);
  }
  if (!(bits & 0x080))
  {
    edges[i  ].first.set( xmin,ymax,zmax);
    edges[i++].second.set(xmin,ymin,zmax);
  }

  if (!(bits & 0x100))
  {
    edges[i  ].first.set( xmin,ymin,zmin);
    edges[i++].second.set(xmin,ymin,zmax);
  }
  if (!(bits & 0x200))
  {
    edges[i  ].first.set( xmax,ymin,zmin);
    edges[i++].second.set(xmax,ymin,zmax);
  }
  if (!(bits & 0x400))
  {
    edges[i  ].first.set( xmax,ymax,zmin);
    edges[i++].second.set(xmax,ymax,zmax);
  }
  if (!(bits & 0x800))
  {
    edges[i  ].first.set( xmin,ymax,zmin);
    edges[i++].second.set(xmin,ymax,zmax);
  }
  edges.resize(i);

  // Clip the edges by the planes
  //
  for (auto iedge = edges.cbegin(); iedge != edges.cend(); ++iedge)
  {
    G4bool    exist = true;
    G4Point3D p1    = iedge->first;
    G4Point3D p2    = iedge->second;
    for (auto iplane = pPlanes.cbegin(); iplane != pPlanes.cend(); ++iplane)
    {
      // Clip current edge
      G4double d1 = iplane->distance(p1);
      G4double d2 = iplane->distance(p2);
      if (d1 > 0.0)
      {
        if (d2 > 0.0) { exist = false; break; } // go to next edge
        p1 = (p2*d1-p1*d2)/(d1-d2);                   // move p1
      }
      else
      {
        if (d2 > 0.0) { p2 = (p1*d2-p2*d1)/(d2-d1); } // move p2
      }
    }
    // Adjust the extent
    if (exist)
    {
      emin.setX(std::min(std::min(p1.x(),p2.x()),emin.x()));
      emin.setY(std::min(std::min(p1.y(),p2.y()),emin.y()));
      emin.setZ(std::min(std::min(p1.z(),p2.z()),emin.z()));

      emax.setX(std::max(std::max(p1.x(),p2.x()),emax.x()));
      emax.setY(std::max(std::max(p1.y(),p2.y()),emax.y()));
      emax.setZ(std::max(std::max(p1.z(),p2.z()),emax.z()));
    }
  }

  // Copy the extent back
  //
  pExtent.first  = emin;
  pExtent.second = emax;
}
