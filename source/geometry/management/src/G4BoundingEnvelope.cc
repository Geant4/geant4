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
//
// Implementation of G4BoundingEnvelope
//
// Author: evgueni.tcherniaev@cern.ch
//
// 2016.05.25 E.Tcherniaev - initial version
// 
// --------------------------------------------------------------------

#include <cmath>
#include "globals.hh"

#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"

///////////////////////////////////////////////////////////////////////
//
// Constructor from an axis aligned bounding box
//
G4BoundingEnvelope::G4BoundingEnvelope(const G4ThreeVector& pMin,
                                       const G4ThreeVector& pMax,
                                       G4double delta)
{
  SetDelta(delta);
  SetBoundingBox(pMin,pMax);
}

///////////////////////////////////////////////////////////////////////
//
// Constructor from a prism
//
G4BoundingEnvelope::G4BoundingEnvelope(const G4ThreeVectorList& baseA,
                                       const G4ThreeVectorList& baseB,
                                       G4double delta)
{
  SetDelta(delta);
  SetBoundingPrism(baseA,baseB);
}

///////////////////////////////////////////////////////////////////////
//
// Constructor from a pyramid
//
G4BoundingEnvelope::G4BoundingEnvelope(const G4ThreeVector& apex,
                                       const G4ThreeVectorList& base,
                                       G4double delta)
{
  SetDelta(delta);
  SetBoundingPyramid(apex,base);
}

///////////////////////////////////////////////////////////////////////
//
// Constructor from a sequence of polygons
//
G4BoundingEnvelope::G4BoundingEnvelope(
  const std::vector<G4ThreeVectorList*>& polygons,G4double delta)
{
  SetDelta(delta);
  SetBoundingPolygons(polygons);
}

///////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4BoundingEnvelope::G4BoundingEnvelope(const G4BoundingEnvelope& rhs)
  : fDelta(rhs.fDelta)
{
  // Copy data
  G4int nb = rhs.fBases.size();
  fBases.resize(nb);
  for (G4int i=0; i<nb; i++) {
    fBases[i] = new G4Polygon3D(*rhs.fBases[i]);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4BoundingEnvelope&
G4BoundingEnvelope::operator=(const G4BoundingEnvelope& rhs)
{
  // Check assignment to self
  if (this == &rhs)  { return *this; }

  // Copy data
  fDelta = rhs.fDelta;
  CleanPolygons();
  G4int nb = rhs.fBases.size();
  fBases.resize(nb);
  for (G4int i=0; i<nb; i++) {
    fBases[i] = new G4Polygon3D(*rhs.fBases[i]);
  }

  return *this;
}

///////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4BoundingEnvelope::~G4BoundingEnvelope()
{
  CleanPolygons();
  fBases.resize(0);
}

///////////////////////////////////////////////////////////////////////
//
// Set the extension
//
void G4BoundingEnvelope::SetDelta(G4double delta)
{
  fDelta = std::abs(delta);
}

///////////////////////////////////////////////////////////////////////
//
// Set axis aligned bounding box
//
void
G4BoundingEnvelope::SetBoundingBox(const G4ThreeVector& pMin,
                                   const G4ThreeVector& pMax)
{
  // Check parameters
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Badly defined bounding box (min >= max)!"
            << "\npMin = " << pMin 
            << "\npMax = " << pMax;
    G4Exception("G4BoundingEnvelope::SetBoundingBox()",
                "GeomMgt0001", FatalException, message);
  }

  CleanPolygons();
  fBases.resize(2);

  // Set 1st base
  fBases[0] = new G4Polygon3D(4);
  (*fBases[0])[0] = G4Point3D(pMin.x(),pMin.y(),pMin.z());
  (*fBases[0])[1] = G4Point3D(pMax.x(),pMin.y(),pMin.z());
  (*fBases[0])[2] = G4Point3D(pMax.x(),pMax.y(),pMin.z());
  (*fBases[0])[3] = G4Point3D(pMin.x(),pMax.y(),pMin.z());

  // Set 2nd base
  fBases[1] = new G4Polygon3D(4);
  (*fBases[1])[0] = G4Point3D(pMin.x(),pMin.y(),pMax.z());
  (*fBases[1])[1] = G4Point3D(pMax.x(),pMin.y(),pMax.z());
  (*fBases[1])[2] = G4Point3D(pMax.x(),pMax.y(),pMax.z());
  (*fBases[1])[3] = G4Point3D(pMin.x(),pMax.y(),pMax.z());
}

///////////////////////////////////////////////////////////////////////
//
// Set bounding prism
//
void
G4BoundingEnvelope::SetBoundingPrism(const G4ThreeVectorList& baseA,
                                     const G4ThreeVectorList& baseB)
{
  G4int na = baseA.size();
  G4int nb = baseB.size();
  if (na < 3 || nb < 3 || na != nb)
  {
    std::ostringstream message;
    message << "Badly defined bases of the bounding prism!"
            << "\nNumber of vertices in 1st base: " << na
            << "\nNumber of vertices in 2nd base: " << nb;
    G4Exception("G4BoundingEnvelope::SetBoundingPrism()",
                "GeomMgt0001", FatalException, message);
  }

  CleanPolygons();
  fBases.resize(2);

  // Set 1st base
  fBases[0] = new G4Polygon3D(na);
  for (G4int i=0; i<na; i++) (*fBases[0])[i] = baseA[i];

  // Set 2nd base
  fBases[1] = new G4Polygon3D(nb);
  for (G4int i=0; i<nb; i++) (*fBases[1])[i] = baseB[i];
}

///////////////////////////////////////////////////////////////////////
//
// Set bounding pyramid
//
void
G4BoundingEnvelope::SetBoundingPyramid(const G4ThreeVector& apex,
                                       const G4ThreeVectorList& base)
{
  // Check parameters
  G4int np = base.size();
  if (np < 3)
  {
    std::ostringstream message;
    message << "Badly defined base of the bounding pyramid!"
            << "\nNumber of vertices in the base: " << np;
    G4Exception("G4BoundingEnvelope::SetBoundingPyramid()", 
                "GeomMgt0001", FatalException, message);
  }

  CleanPolygons();
  fBases.resize(2);

  // Set apex
  fBases[0] = new G4Polygon3D(1);
  (*fBases[0])[0] = apex;

  // Set base
  fBases[1] = new G4Polygon3D(np);
  for (G4int i=0; i<np; i++) (*fBases[1])[i] = base[i];
}

///////////////////////////////////////////////////////////////////////
//
// Set bounding sequence of polygons.
// Firsf and last polygons may consist of a single vertex
//
void G4BoundingEnvelope::SetBoundingPolygons(
       const std::vector<G4ThreeVectorList*>& polygons)
{
  // Check parameters
  G4int nbases = polygons.size();
  if (nbases < 2)
  {
    std::ostringstream message;
    message << "Wrong number of polygons in the sequence: " << nbases
            << "\nShould be at least two!";
    G4Exception("G4BoundingEnvelope::SetBoundingPolygons()", 
                "GeomMgt0001", FatalException, message);
    return;
  }

  G4int nsize  = std::max(polygons[0]->size(),polygons[1]->size());
  if (nsize < 3) {
    std::ostringstream message;
    message << "Badly constructed polygons!"
            << "\nNumber of polygons: " << nbases
            << "\nPolygon #0 size: " << polygons[0]->size()
            << "\nPolygon #1 size: " << polygons[1]->size()
            << "\n...";
    G4Exception("G4BoundingEnvelope::SetBoundingPolygons()", 
                "GeomMgt0001", FatalException, message);
    return;
  }

  for (G4int k=0; k<nbases; k++) {
    G4int np = polygons[k]->size();
    if (np == nsize)          continue;
    if (np == 1 && k==0)      continue;
    if (np == 1 && k==nbases) continue;
    std::ostringstream message;
    message << "Badly constructed polygons!"
            << "\nNumber of polygons: " << nbases
            << "\nPolygon #" << k << " size: " << np
            << "\nexpected size: " << nsize;
    G4Exception("G4BoundingEnvelope::SetBoundingPolygons()", 
                "GeomMgt0001", FatalException, message);
    return;
  }  

  // Copy polygons
  CleanPolygons();
  fBases.resize(nbases);
  for (G4int k=0; k<nbases; k++) {
    G4int np = polygons[k]->size();
    fBases[k] = new G4Polygon3D(np); 
    for (G4int i=0; i<np; i++) (*fBases[k])[i] = (*polygons[k])[i];
  }
}

///////////////////////////////////////////////////////////////////////
//
// Free memory allocated for polygons
//
void G4BoundingEnvelope::CleanPolygons()
{
  G4int nb = fBases.size();
  for (G4int i=0; i<nb; i++) {
    delete fBases[i]; fBases[i] = 0;
  }
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
  // Create adjusted G4VoxelLimits box. New limits are extended by
  // fDelta multiplied by max scale factor of the transformation.
  //
  G4Scale3D scale3D; G4Rotate3D rotate3D; G4Translate3D translate3D;
  pTransform3D.getDecomposition(scale3D, rotate3D, translate3D);
  G4double scale = std::max(std::max(std::abs(scale3D.xx()),
                                     std::abs(scale3D.yy())),
                                     std::abs(scale3D.zz()));
  G4double delta = (scale > 1.) ? fDelta*scale : fDelta;
  G4VoxelLimits limits = GetAdjustedVoxelLimits(pVoxelLimits, delta);

  // Main loop along the set of prisms
  //
  G4Segment3D extent(G4Point3D( kInfinity, kInfinity, kInfinity),
                     G4Point3D(-kInfinity,-kInfinity,-kInfinity));

  G4int nbases = fBases.size();
  for (G4int k=0; k<nbases-1; k++)
  {
    // Transform vertices of and find bounding box of current prism
    G4Polygon3D baseA, baseB;
    G4Segment3D prismAABB(G4Point3D( kInfinity, kInfinity, kInfinity),
                          G4Point3D(-kInfinity,-kInfinity,-kInfinity));

    TransformVertices(pTransform3D, *fBases[k]  , baseA, prismAABB);
    TransformVertices(pTransform3D, *fBases[k+1], baseB, prismAABB);

    // Check that bounding box of the prism intersect the voxel limits
    if (prismAABB.first.x()  > limits.GetMaxXExtent()) continue;
    if (prismAABB.first.y()  > limits.GetMaxYExtent()) continue;
    if (prismAABB.first.z()  > limits.GetMaxZExtent()) continue;
    if (prismAABB.second.x() < limits.GetMinXExtent()) continue;
    if (prismAABB.second.y() < limits.GetMinYExtent()) continue;
    if (prismAABB.second.z() < limits.GetMinZExtent()) continue;

    // Clip edges of the prism by adjusted G4VoxelLimits box
    std::vector<G4Segment3D> vecEdges; 
    CreateListOfEdges(baseA, baseB, vecEdges);
    if (ClipEdgesByVoxelLimits(vecEdges, limits, extent)) continue;

    // Some edges of the prism are completely outside of the voxel
    // limits, clip edges of adjusted G4VoxelLimits box by the prism 
    std::vector<G4Plane3D> vecPlanes; 
    CreateListOfPlanes(baseA, baseB, vecPlanes);
    ClipVoxelLimitsByPlanes(limits, vecPlanes, prismAABB, extent);
  }

  // Final adjustment of the extent
  // 
  G4double emin=kInfinity, emax=kInfinity;
  if (pAxis == kXAxis) { emin = extent.first.x(); emax = extent.second.x(); }
  if (pAxis == kYAxis) { emin = extent.first.y(); emax = extent.second.y(); }
  if (pAxis == kZAxis) { emin = extent.first.z(); emax = extent.second.z(); }

  G4bool exist = false;
  if (emin <= emax) {
    exist = true;
    // Add the extension to the endpoints
    if (emin > limits.GetMinExtent(pAxis)) emin -= delta;
    if (emax < limits.GetMaxExtent(pAxis)) emax += delta;

    G4double kCarTolerance = 
      G4GeometryTolerance::GetInstance()->GetSurfaceTolerance(); 

    // Clip by original voxel limits, if required
    if (emin <= pVoxelLimits.GetMinExtent(pAxis)) {
      pMin = pVoxelLimits.GetMinExtent(pAxis) - kCarTolerance;
    } else {
      pMin = emin;
    }
    if (emax >= pVoxelLimits.GetMaxExtent(pAxis)) {
      pMax = pVoxelLimits.GetMaxExtent(pAxis) + kCarTolerance;
    } else {
      pMax = emax;
    }
    exist = true;
  } else {
    exist = false;
    pMin =  kInfinity;
    pMax = -kInfinity;
  }

  return exist;
}

///////////////////////////////////////////////////////////////////////
//
// Create adjusted voxel limits
//
G4VoxelLimits
G4BoundingEnvelope::GetAdjustedVoxelLimits(const G4VoxelLimits& pVoxelLimits,
                                           G4double pDelta) const
{
  EAxis axis[] = { kXAxis,kYAxis,kZAxis };
  G4VoxelLimits limits; // default is unlimited 
  for (G4int i=0; i<3; i++) {
    if (pVoxelLimits.IsLimited(axis[i])) {
      G4double emin = pVoxelLimits.GetMinExtent(axis[i]) - pDelta;
      G4double emax = pVoxelLimits.GetMaxExtent(axis[i]) + pDelta;
      limits.AddLimit(axis[i], emin, emax);
   }
  }
  return limits;
}

///////////////////////////////////////////////////////////////////////
//
// Transform vertices of a polygon and update the bounding box
//
void
G4BoundingEnvelope::TransformVertices(const G4Transform3D& pTransform3D,
                                      const G4Polygon3D& polyA,
                                            G4Polygon3D& polyB,
                                            G4Segment3D& pAABB) const
{
  G4double xmin = pAABB.first.x();
  G4double ymin = pAABB.first.y();
  G4double zmin = pAABB.first.z();
  G4double xmax = pAABB.second.x();
  G4double ymax = pAABB.second.y();
  G4double zmax = pAABB.second.z();

  G4int np = polyA.size();
  polyB.resize(np);
  for (G4int i=0; i<np; i++) {
    polyB[i] = pTransform3D*polyA[i];
    xmin = std::min(xmin,polyB[i].x());
    ymin = std::min(ymin,polyB[i].y());
    zmin = std::min(zmin,polyB[i].z());
    xmax = std::max(xmax,polyB[i].x());
    ymax = std::max(ymax,polyB[i].y());
    zmax = std::max(zmax,polyB[i].z());
  }

  pAABB.first.set( xmin,ymin,zmin);
  pAABB.second.set(xmax,ymax,zmax);
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
  G4int na = baseA.size();
  G4int nb = baseB.size();
  pEdges.resize(0);
  if (na == nb) {
    G4int k = na - 1;
    for (G4int i=0; i<na; i++) {
      pEdges.push_back(G4Segment3D(baseA[i],baseB[i]));
      pEdges.push_back(G4Segment3D(baseA[i],baseA[k]));
      pEdges.push_back(G4Segment3D(baseB[i],baseB[k]));
      k = i;
    }
  } else if (nb == 1) {
    G4int k = na - 1;
    for (G4int i=0; i<na; i++) {
      pEdges.push_back(G4Segment3D(baseA[i],baseA[k]));
      pEdges.push_back(G4Segment3D(baseA[i],baseB[0]));
      k = i;
    }
  } else if (na == 1) {
    G4int k = nb - 1;
    for (G4int i=0; i<nb; i++) {
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
  G4int na = baseA.size();
  G4int nb = baseB.size();
  G4Point3D pa(0.,0.,0.), pb(0.,0.,0.), p0;
  for (G4int i=0; i<na; i++) pa += baseA[i];
  for (G4int i=0; i<nb; i++) pb += baseB[i];
  pa /= na; pb /= nb; p0 = (pa+pb)/2.;

  // Create list of planes
  //
  pPlanes.resize(0);
  if (na == nb) {
    G4int k = na - 1;
    for (G4int i=0; i<na; i++) {
      pPlanes.push_back(G4Plane3D(baseA[i],baseA[k],baseB[k]));
      k = i;
    }
  pPlanes.push_back(G4Plane3D(baseA[1],baseA[0],pa));
  pPlanes.push_back(G4Plane3D(baseB[0],baseB[1],pb));
  } else if (nb == 1) {
    G4int k = na - 1;
    for (G4int i=0; i<na; i++) {
      pPlanes.push_back(G4Plane3D(baseA[i],baseA[k],baseB[0]));
      k = i;
    }
    pPlanes.push_back(G4Plane3D(baseA[2],baseA[1],baseA[0]));
  } else if (na == 1) {
    G4int k = nb - 1;
    for (G4int i=0; i<nb; i++) {
      pPlanes.push_back(G4Plane3D(baseB[k],baseB[i],baseA[0]));
      k = i;
    }
    pPlanes.push_back(G4Plane3D(baseB[0],baseB[1],baseB[2]));
  }

  // Ensure that normals of the planes point to outside
  //
  G4int nplanes = pPlanes.size();
  for (G4int i=0; i<nplanes; i++) {
    pPlanes[i].normalize();
    if (pPlanes[i].distance(p0) > 0) {
      pPlanes[i] = G4Plane3D(-pPlanes[i].a(),-pPlanes[i].b(),
                             -pPlanes[i].c(),-pPlanes[i].d());
    } 
  }
}

///////////////////////////////////////////////////////////////////////
//
// Clip edges of a prism by G4VoxelLimits box  
//
G4bool
G4BoundingEnvelope::ClipEdgesByVoxelLimits(const std::vector<G4Segment3D>& pEdges,
                                           const G4VoxelLimits& pBox,
                                           G4Segment3D& pExtent) const
{
  G4bool    done = true;
  G4Point3D emin = pExtent.first;
  G4Point3D emax = pExtent.second;

  G4int nedges = pEdges.size();
  for (G4int k=0; k<nedges; k++)
  {
    G4double  d1, d2;
    G4Point3D p1 = pEdges[k].first;
    G4Point3D p2 = pEdges[k].second;

    // Clip current edge by X min
    d1 = pBox.GetMinXExtent() - p1.x(); 
    d2 = pBox.GetMinXExtent() - p2.x(); 
    if (d1 > 0.0) {
      if (d2 > 0.0) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);                   // move p1
    } else {
      if (d2 > 0.0) { p2 = (p1*d2-p2*d1)/(d2-d1); } // move p2
    }

    // Clip current edge by X max
    d1 = p1.x() - pBox.GetMaxXExtent(); 
    d2 = p2.x() - pBox.GetMaxXExtent(); 
    if (d1 > 0.) {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    } else {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Y min
    d1 = pBox.GetMinYExtent() - p1.y(); 
    d2 = pBox.GetMinYExtent() - p2.y(); 
    if (d1 > 0.) {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    } else {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Y max
    d1 = p1.y() - pBox.GetMaxYExtent(); 
    d2 = p2.y() - pBox.GetMaxYExtent(); 
    if (d1 > 0.) {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    } else {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Z min
    d1 = pBox.GetMinZExtent() - p1.z(); 
    d2 = pBox.GetMinZExtent() - p2.z(); 
    if (d1 > 0.) {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    } else {
      if (d2 > 0.) { p2 = (p1*d2-p2*d1)/(d2-d1); }
    }

    // Clip current edge by Z max
    d1 = p1.z() - pBox.GetMaxZExtent(); 
    d2 = p2.z() - pBox.GetMaxZExtent(); 
    if (d1 > 0.) {
      if (d2 > 0.) { done = false; continue; } // go to next edge
      p1 = (p2*d1-p1*d2)/(d1-d2);
    } else {
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
// Clip G4VoxelLimits by set of planes bounding a convex prism
//
void
G4BoundingEnvelope::ClipVoxelLimitsByPlanes(const G4VoxelLimits& pBox,
                                            const std::vector<G4Plane3D>& pPlanes,
                                            const G4Segment3D& pAABB,
                                                  G4Segment3D& pExtent) const
{
  G4Point3D emin = pExtent.first;
  G4Point3D emax = pExtent.second;

  // Create 12 edges of the voxel limits box, reduce them where
  // appropriate to avoid calculations with big numbers (kInfinity)
  //
  G4double xmin = pBox.GetMinXExtent(), xmax = pBox.GetMaxXExtent();
  G4double ymin = pBox.GetMinYExtent(), ymax = pBox.GetMaxYExtent();
  G4double zmin = pBox.GetMinZExtent(), zmax = pBox.GetMaxZExtent();
  if( xmin < 2.*pAABB.first.x() && xmax > 2.*pAABB.second.x())
    { xmin = 2.*pAABB.first.x();   xmax = 2.*pAABB.second.x(); }
  if( ymin < 2.*pAABB.first.y() && ymax > 2.*pAABB.second.y())
    { ymin = 2.*pAABB.first.y();   ymax = 2.*pAABB.second.y(); }
  if( zmin < 2.*pAABB.first.z() && zmax > 2.*pAABB.second.z())
    { zmin = 2.*pAABB.first.z();   zmax = 2.*pAABB.second.z(); }

  std::vector<G4Segment3D> edges(12);
  edges[0].first.set(xmin,ymin,zmin); edges[0].second.set(xmax,ymin,zmin);
  edges[1].first = edges[0].second;   edges[1].second.set(xmax,ymax,zmin);
  edges[2].first = edges[1].second;   edges[2].second.set(xmin,ymax,zmin);
  edges[3].first = edges[2].second;   edges[3].second = edges[0].first;

  edges[4].first.set(xmin,ymin,zmax); edges[4].second.set(xmax,ymin,zmax);
  edges[5].first = edges[4].second;   edges[5].second.set(xmax,ymax,zmax);
  edges[6].first = edges[5].second;   edges[6].second.set(xmin,ymax,zmax);
  edges[7].first = edges[6].second;   edges[7].second = edges[4].first;

  edges[ 8].first = edges[0].first;   edges[ 8].second = edges[4].first;
  edges[ 9].first = edges[1].first;   edges[ 9].second = edges[5].first;
  edges[10].first = edges[2].first;   edges[10].second = edges[6].first;
  edges[11].first = edges[3].first;   edges[11].second = edges[7].first;

  // Clip the edges by the planes
  //
  G4int nedges  = edges.size();
  G4int nplanes = pPlanes.size();
  for (G4int k=0; k<nedges; k++)
  {
    G4Point3D p1 = edges[k].first;
    G4Point3D p2 = edges[k].second;
    G4bool exist = true;
    for (G4int i=0; i<nplanes; i++) {
      // Clip current edge
      G4double d1 = pPlanes[i].distance(p1); 
      G4double d2 = pPlanes[i].distance(p2); 
      if (d1 > 0.0) {
        if (d2 > 0.0) { exist = false; break; } // go to next edge
        p1 = (p2*d1-p1*d2)/(d1-d2);                   // move p1
      } else {
        if (d2 > 0.0) { p2 = (p1*d2-p2*d1)/(d2-d1); } // move p2
      }
    }
    // Adjust the extent
    if (exist) {
      emin.setX(std::min(std::min(p1.x(),p2.x()),emin.x())); 
      emin.setY(std::min(std::min(p1.y(),p2.y()),emin.y())); 
      emin.setZ(std::min(std::min(p1.z(),p2.z()),emin.z())); 

      emax.setX(std::max(std::max(p1.x(),p2.x()),emax.x())); 
      emax.setY(std::max(std::max(p1.y(),p2.y()),emax.y())); 
      emax.setZ(std::max(std::max(p1.z(),p2.z()),emax.z())); 
    }
  }

  // Copy the extent back
  pExtent.first  = emin;
  pExtent.second = emax;
}
