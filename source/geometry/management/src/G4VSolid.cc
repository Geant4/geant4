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
// G4VSolid implementation for solid base class
//
// 10.10.18 E.Tcherniaev, more robust EstimateSurfaceArea() based on distance
// 30.06.95 P.Kent, Created.
// --------------------------------------------------------------------

#include "G4VSolid.hh"
#include "G4SolidStore.hh"
#include "globals.hh"
#include "G4QuickRand.hh"
#include "G4GeometryTolerance.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VisExtent.hh"

//////////////////////////////////////////////////////////////////////////
//
// Streaming operator dumping solid contents

std::ostream& operator<< ( std::ostream& os, const G4VSolid& e )
{
    return e.StreamInfo(os);
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor
//  - Copies name
//  - Add ourselves to solid Store

G4VSolid::G4VSolid(const G4String& name)
  : fshapeName(name)
{
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

    // Register to store
    //
    G4SolidStore::GetInstance()->Register(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//

G4VSolid::G4VSolid(const G4VSolid& rhs)
  : kCarTolerance(rhs.kCarTolerance), fshapeName(rhs.fshapeName)
{
    // Register to store
    //
    G4SolidStore::GetInstance()->Register(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4VSolid::G4VSolid( __void__& )
  : fshapeName("")
{
    // Register to store
    //
    G4SolidStore::GetInstance()->Register(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor (virtual)
// - Remove ourselves from solid Store

G4VSolid::~G4VSolid()
{
    G4SolidStore::GetInstance()->DeRegister(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4VSolid& G4VSolid::operator = (const G4VSolid& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy data
   //
   kCarTolerance = rhs.kCarTolerance;
   fshapeName = rhs.fshapeName;

   return *this;
}



//////////////////////////////////////////////////////////////////////////
//
// Set solid name and notify store of the change

void G4VSolid::SetName(const G4String& name)
{
  fshapeName = name;
  G4SolidStore::GetInstance()->SetMapValid(false);
}

//////////////////////////////////////////////////////////////////////////
//
// Throw exception if ComputeDimensions called for illegal derived class

void G4VSolid::ComputeDimensions(G4VPVParameterisation*,
                                 const G4int,
                                 const G4VPhysicalVolume*)
{
    std::ostringstream message;
    message << "Illegal call to G4VSolid::ComputeDimensions()" << G4endl
            << "Method not overloaded by derived class !";
    G4Exception("G4VSolid::ComputeDimensions()", "GeomMgt0003",
                FatalException, message);
}

//////////////////////////////////////////////////////////////////////////
//
// Throw exception (warning) for solids not implementing the method

G4ThreeVector G4VSolid::GetPointOnSurface() const
{
    std::ostringstream message;
    message << "Not implemented for solid: "
            << GetEntityType() << " !" << G4endl
            << "Returning origin.";
    G4Exception("G4VSolid::GetPointOnSurface()", "GeomMgt1001",
                JustWarning, message);
    return G4ThreeVector(0,0,0);
}

//////////////////////////////////////////////////////////////////////////
//
// Dummy implementations ...

const G4VSolid* G4VSolid::GetConstituentSolid(G4int) const
{ return nullptr; }

G4VSolid* G4VSolid::GetConstituentSolid(G4int)
{ return nullptr; }

const G4DisplacedSolid* G4VSolid::GetDisplacedSolidPtr() const
{ return nullptr; }

G4DisplacedSolid* G4VSolid::GetDisplacedSolidPtr()
{ return nullptr; }

////////////////////////////////////////////////////////////////
//
// Returns an estimation of the solid volume in internal units.
// The number of statistics and error accuracy is fixed.
// This method may be overloaded by derived classes to compute the
// exact geometrical quantity for solids where this is possible.
// or anyway to cache the computed value.
// This implementation does NOT cache the computed value.

G4double G4VSolid::GetCubicVolume()
{
  G4int cubVolStatistics = 1000000;
  G4double cubVolEpsilon = 0.001;
  return EstimateCubicVolume(cubVolStatistics, cubVolEpsilon);
}

////////////////////////////////////////////////////////////////
//
// Calculate cubic volume based on Inside() method.
// Accuracy is limited by the second argument or the statistics
// expressed by the first argument.
// Implementation is courtesy of Vasiliki Despoina Mitsou,
// University of Athens.

G4double G4VSolid::EstimateCubicVolume(G4int nStat, G4double epsilon) const
{
  G4int iInside=0;
  G4double px,py,pz,minX,maxX,minY,maxY,minZ,maxZ,volume,halfepsilon;
  G4ThreeVector p;
  EInside in;

  // values needed for CalculateExtent signature

  G4VoxelLimits limit;                // Unlimited
  G4AffineTransform origin;

  // min max extents of pSolid along X,Y,Z

  CalculateExtent(kXAxis,limit,origin,minX,maxX);
  CalculateExtent(kYAxis,limit,origin,minY,maxY);
  CalculateExtent(kZAxis,limit,origin,minZ,maxZ);

  // limits

  if(nStat < 100)    nStat   = 100;
  if(epsilon > 0.01) epsilon = 0.01;
  halfepsilon = 0.5*epsilon;

  for(auto i = 0; i < nStat; ++i )
  {
    px = minX-halfepsilon+(maxX-minX+epsilon)*G4QuickRand();
    py = minY-halfepsilon+(maxY-minY+epsilon)*G4QuickRand();
    pz = minZ-halfepsilon+(maxZ-minZ+epsilon)*G4QuickRand();
    p  = G4ThreeVector(px,py,pz);
    in = Inside(p);
    if(in != kOutside) ++iInside;
  }
  volume = (maxX-minX+epsilon)*(maxY-minY+epsilon)
         * (maxZ-minZ+epsilon)*iInside/nStat;
  return volume;
}

////////////////////////////////////////////////////////////////
//
// Returns an estimation of the solid surface area in internal units.
// The number of statistics and error accuracy is fixed.
// This method may be overloaded by derived classes to compute the
// exact geometrical quantity for solids where this is possible.
// or anyway to cache the computed value.
// This implementation does NOT cache the computed value.

G4double G4VSolid::GetSurfaceArea()
{
  G4int stat = 1000000;
  G4double ell = -1.;
  return EstimateSurfaceArea(stat,ell);
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate surface area by estimating volume of a thin shell
// surrounding the surface using Monte-Carlo method.
// Input parameters:
//    nstat - statistics (number of random points)
//    eps   - shell thinkness

G4double G4VSolid::EstimateSurfaceArea(G4int nstat, G4double ell) const
{
  static const G4double s2 = 1./std::sqrt(2.);
  static const G4double s3 = 1./std::sqrt(3.);
  static const G4ThreeVector directions[64] =
  {
    G4ThreeVector(  0,  0,  0), G4ThreeVector( -1,  0,  0), // (  ,  ,  ) ( -,  ,  )
    G4ThreeVector(  1,  0,  0), G4ThreeVector( -1,  0,  0), // ( +,  ,  ) (-+,  ,  )
    G4ThreeVector(  0, -1,  0), G4ThreeVector(-s2,-s2,  0), // (  , -,  ) ( -, -,  )
    G4ThreeVector( s2, -s2, 0), G4ThreeVector(  0, -1,  0), // ( +, -,  ) (-+, -,  )

    G4ThreeVector(  0,  1,  0), G4ThreeVector( -s2, s2, 0), // (  , +,  ) ( -, +,  )
    G4ThreeVector( s2, s2,  0), G4ThreeVector(  0,  1,  0), // ( +, +,  ) (-+, +,  )
    G4ThreeVector(  0, -1,  0), G4ThreeVector( -1,  0,  0), // (  ,-+,  ) ( -,-+,  )
    G4ThreeVector(  1,  0,  0), G4ThreeVector( -1,  0,  0), // ( +,-+,  ) (-+,-+,  )

    G4ThreeVector(  0,  0, -1), G4ThreeVector(-s2,  0,-s2), // (  ,  , -) ( -,  , -)
    G4ThreeVector( s2,  0,-s2), G4ThreeVector(  0,  0, -1), // ( +,  , -) (-+,  , -)
    G4ThreeVector(  0,-s2,-s2), G4ThreeVector(-s3,-s3,-s3), // (  , -, -) ( -, -, -)
    G4ThreeVector( s3,-s3,-s3), G4ThreeVector(  0,-s2,-s2), // ( +, -, -) (-+, -, -)

    G4ThreeVector(  0, s2,-s2), G4ThreeVector(-s3, s3,-s3), // (  , +, -) ( -, +, -)
    G4ThreeVector( s3, s3,-s3), G4ThreeVector(  0, s2,-s2), // ( +, +, -) (-+, +, -)
    G4ThreeVector(  0,  0, -1), G4ThreeVector(-s2,  0,-s2), // (  ,-+, -) ( -,-+, -)
    G4ThreeVector( s2,  0,-s2), G4ThreeVector(  0,  0, -1), // ( +,-+, -) (-+,-+, -)

    G4ThreeVector(  0,  0,  1), G4ThreeVector(-s2,  0, s2), // (  ,  , +) ( -,  , +)
    G4ThreeVector( s2,  0, s2), G4ThreeVector(  0,  0,  1), // ( +,  , +) (-+,  , +)
    G4ThreeVector(  0,-s2, s2), G4ThreeVector(-s3,-s3, s3), // (  , -, +) ( -, -, +)
    G4ThreeVector( s3,-s3, s3), G4ThreeVector(  0,-s2, s2), // ( +, -, +) (-+, -, +)

    G4ThreeVector(  0, s2, s2), G4ThreeVector(-s3, s3, s3), // (  , +, +) ( -, +, +)
    G4ThreeVector( s3, s3, s3), G4ThreeVector(  0, s2, s2), // ( +, +, +) (-+, +, +)
    G4ThreeVector(  0,  0,  1), G4ThreeVector(-s2,  0, s2), // (  ,-+, +) ( -,-+, +)
    G4ThreeVector( s2,  0, s2), G4ThreeVector(  0,  0,  1), // ( +,-+, +) (-+,-+, +)

    G4ThreeVector(  0,  0, -1), G4ThreeVector( -1,  0,  0), // (  ,  ,-+) ( -,  ,-+)
    G4ThreeVector(  1,  0,  0), G4ThreeVector( -1,  0,  0), // ( +,  ,-+) (-+,  ,-+)
    G4ThreeVector(  0, -1,  0), G4ThreeVector(-s2,-s2,  0), // (  , -,-+) ( -, -,-+)
    G4ThreeVector( s2, -s2, 0), G4ThreeVector(  0, -1,  0), // ( +, -,-+) (-+, -,-+)

    G4ThreeVector(  0,  1,  0), G4ThreeVector( -s2, s2, 0), // (  , +,-+) ( -, +,-+)
    G4ThreeVector( s2, s2,  0), G4ThreeVector(  0,  1,  0), // ( +, +,-+) (-+, +,-+)
    G4ThreeVector(  0, -1,  0), G4ThreeVector( -1,  0,  0), // (  ,-+,-+) ( -,-+,-+)
    G4ThreeVector(  1,  0,  0), G4ThreeVector( -1,  0,  0), // ( +,-+,-+) (-+,-+,-+)
  };

  G4ThreeVector bmin, bmax;
  BoundingLimits(bmin, bmax);

  G4double dX = bmax.x() - bmin.x();
  G4double dY = bmax.y() - bmin.y();
  G4double dZ = bmax.z() - bmin.z();

  // Define statistics and shell thickness
  //
  G4int npoints = (nstat < 1000) ? 1000 : nstat;
  G4double coeff = 0.5 / std::cbrt(G4double(npoints));
  G4double eps = (ell > 0) ? ell : coeff * std::min(std::min(dX, dY), dZ);
  G4double del = 1.8 * eps; // shold be more than sqrt(3.)

  G4double minX = bmin.x() - eps;
  G4double minY = bmin.y() - eps;
  G4double minZ = bmin.z() - eps;

  G4double dd = 2. * eps;
  dX += dd;
  dY += dd;
  dZ += dd;

  // Calculate surface area
  //
  G4int icount = 0;
  for(auto i = 0; i < npoints; ++i)
  {
    G4double px = minX + dX*G4QuickRand();
    G4double py = minY + dY*G4QuickRand();
    G4double pz = minZ + dZ*G4QuickRand();
    G4ThreeVector p  = G4ThreeVector(px, py, pz);
    EInside in = Inside(p);
    G4double dist = 0;
    if (in == kInside)
    {
      if (DistanceToOut(p) >= eps) continue;
      G4int icase = 0;
      if (Inside(G4ThreeVector(px-del, py, pz)) != kInside) icase += 1;
      if (Inside(G4ThreeVector(px+del, py, pz)) != kInside) icase += 2;
      if (Inside(G4ThreeVector(px, py-del, pz)) != kInside) icase += 4;
      if (Inside(G4ThreeVector(px, py+del, pz)) != kInside) icase += 8;
      if (Inside(G4ThreeVector(px, py, pz-del)) != kInside) icase += 16;
      if (Inside(G4ThreeVector(px, py, pz+del)) != kInside) icase += 32;
      if (icase == 0) continue;
      G4ThreeVector v = directions[icase];
      dist = DistanceToOut(p, v);
      G4ThreeVector n = SurfaceNormal(p + v*dist);
      dist *= v.dot(n);
    }
    else if (in == kOutside)
    {
      if (DistanceToIn(p) >= eps) continue;
      G4int icase = 0;
      if (Inside(G4ThreeVector(px-del, py, pz)) != kOutside) icase += 1;
      if (Inside(G4ThreeVector(px+del, py, pz)) != kOutside) icase += 2;
      if (Inside(G4ThreeVector(px, py-del, pz)) != kOutside) icase += 4;
      if (Inside(G4ThreeVector(px, py+del, pz)) != kOutside) icase += 8;
      if (Inside(G4ThreeVector(px, py, pz-del)) != kOutside) icase += 16;
      if (Inside(G4ThreeVector(px, py, pz+del)) != kOutside) icase += 32;
      if (icase == 0) continue;
      G4ThreeVector v = directions[icase];
      dist = DistanceToIn(p, v);
      if (dist == kInfinity) continue;
      G4ThreeVector n = SurfaceNormal(p + v*dist);
      dist *= -(v.dot(n));
    }
    if (dist < eps) ++icount;
  }
  return dX*dY*dZ*icount/npoints/dd;
}

///////////////////////////////////////////////////////////////////////////
//
// Returns a pointer of a dynamically allocated copy of the solid.
// Returns NULL pointer with warning in case the concrete solid does not
// implement this method. The caller has responsibility for ownership.
//

G4VSolid* G4VSolid::Clone() const
{
  std::ostringstream message;
  message << "Clone() method not implemented for type: "
          << GetEntityType() << "!" << G4endl
          << "Returning NULL pointer!";
  G4Exception("G4VSolid::Clone()", "GeomMgt1001", JustWarning, message);
  return nullptr;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the polygon described
// by the vertices: pSectionIndex->pSectionIndex+1->
//                   pSectionIndex+2->pSectionIndex+3->pSectionIndex
// in the List pVertices
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices
//

void G4VSolid::ClipCrossSection(       G4ThreeVectorList* pVertices,
                                 const G4int pSectionIndex,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const EAxis pAxis,
                                       G4double& pMin, G4double& pMax) const
{

  G4ThreeVectorList polygon;
  polygon.reserve(4);
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  return;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the polygons
// joining the CrossSections at pSectionIndex->pSectionIndex+3 and
//                              pSectionIndex+4->pSectionIndex7
//
// in the List pVertices, within the boundaries of the voxel limits pVoxelLimit
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices

void G4VSolid::ClipBetweenSections(      G4ThreeVectorList* pVertices,
                                   const G4int pSectionIndex,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const EAxis pAxis,
                                         G4double& pMin, G4double& pMax) const
{
  G4ThreeVectorList polygon;
  polygon.reserve(4);
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+3]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex]);
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  return;
}


///////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the convex polygon pPolygon
// along the axis pAxis, within the limits pVoxelLimit
//

void
G4VSolid::CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const EAxis pAxis,
                                        G4double& pMin,
                                        G4double& pMax) const
{
  G4int noLeft,i;
  G4double component;

  ClipPolygon(pPolygon,pVoxelLimit,pAxis);
  noLeft = (G4int)pPolygon.size();

  if ( noLeft )
  {
    for (i=0; i<noLeft; ++i)
    {
      component = pPolygon[i].operator()(pAxis);

      if (component < pMin)
      {
        pMin = component;
      }
      if (component > pMax)
      {
        pMax = component;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Clip the convex polygon described by the vertices at
// pSectionIndex ->pSectionIndex+3 within pVertices to the limits pVoxelLimit
//
// Set pMin to the smallest
//
// Calculate the extent of the polygon along pAxis, when clipped to the
// limits pVoxelLimit. If the polygon exists after clippin, set pMin to
// the polygon's minimum extent along the axis if <pMin, and set pMax to
// the polygon's maximum extent along the axis if >pMax.
//
// The polygon is described by a set of vectors, where each vector represents
// a vertex, so that the polygon is described by the vertex sequence:
//   0th->1st 1st->2nd 2nd->... nth->0th
//
// Modifications to the polygon are made
//
// NOTE: Execessive copying during clipping

void G4VSolid::ClipPolygon(      G4ThreeVectorList& pPolygon,
                           const G4VoxelLimits& pVoxelLimit,
                           const EAxis                        ) const
{
  G4ThreeVectorList outputPolygon;

  if ( pVoxelLimit.IsLimited() )
  {
    if (pVoxelLimit.IsXLimited() ) // && pAxis != kXAxis)
    {
      G4VoxelLimits simpleLimit1;
      simpleLimit1.AddLimit(kXAxis,pVoxelLimit.GetMinXExtent(),kInfinity);
      ClipPolygonToSimpleLimits(pPolygon,outputPolygon,simpleLimit1);

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;

      G4VoxelLimits simpleLimit2;
      simpleLimit2.AddLimit(kXAxis,-kInfinity,pVoxelLimit.GetMaxXExtent());
      ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

      if ( !pPolygon.size() )       return;
      else                          outputPolygon.clear();
    }
    if ( pVoxelLimit.IsYLimited() ) // && pAxis != kYAxis)
    {
      G4VoxelLimits simpleLimit1;
      simpleLimit1.AddLimit(kYAxis,pVoxelLimit.GetMinYExtent(),kInfinity);
      ClipPolygonToSimpleLimits(pPolygon,outputPolygon,simpleLimit1);

      // Must always clear pPolygon - for clip to simpleLimit2 and in case of
      // early exit

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;

      G4VoxelLimits simpleLimit2;
      simpleLimit2.AddLimit(kYAxis,-kInfinity,pVoxelLimit.GetMaxYExtent());
      ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

      if ( !pPolygon.size() )       return;
      else                          outputPolygon.clear();
    }
    if ( pVoxelLimit.IsZLimited() ) // && pAxis != kZAxis)
    {
      G4VoxelLimits simpleLimit1;
      simpleLimit1.AddLimit(kZAxis,pVoxelLimit.GetMinZExtent(),kInfinity);
      ClipPolygonToSimpleLimits(pPolygon,outputPolygon,simpleLimit1);

      // Must always clear pPolygon - for clip to simpleLimit2 and in case of
      // early exit

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;

      G4VoxelLimits simpleLimit2;
      simpleLimit2.AddLimit(kZAxis,-kInfinity,pVoxelLimit.GetMaxZExtent());
      ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

      // Return after final clip - no cleanup
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//
// pVoxelLimits must be only limited along one axis, and either the maximum
// along the axis must be +kInfinity, or the minimum -kInfinity

void
G4VSolid::ClipPolygonToSimpleLimits( G4ThreeVectorList& pPolygon,
                                     G4ThreeVectorList& outputPolygon,
                               const G4VoxelLimits& pVoxelLimit       ) const
{
  G4int i;
  G4int noVertices = (G4int)pPolygon.size();
  G4ThreeVector vEnd,vStart;

  for (i = 0 ; i < noVertices ; ++i )
  {
    vStart = pPolygon[i];
    if ( i == noVertices-1 )    vEnd = pPolygon[0];
    else                        vEnd = pPolygon[i+1];

    if ( pVoxelLimit.Inside(vStart) )
    {
      if (pVoxelLimit.Inside(vEnd))
      {
        // vStart and vEnd inside -> output end point
        //
        outputPolygon.push_back(vEnd);
      }
      else
      {
        // vStart inside, vEnd outside -> output crossing point
        //
        pVoxelLimit.ClipToLimits(vStart,vEnd);
        outputPolygon.push_back(vEnd);
      }
    }
    else
    {
      if (pVoxelLimit.Inside(vEnd))
      {
        // vStart outside, vEnd inside -> output inside section
        //
        pVoxelLimit.ClipToLimits(vStart,vEnd);
        outputPolygon.push_back(vStart);
        outputPolygon.push_back(vEnd);
      }
      else  // Both point outside -> no output
      {
        // outputPolygon.push_back(vStart);
        // outputPolygon.push_back(vEnd);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Throw exception (warning) for solids not implementing the method

void G4VSolid::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  std::ostringstream message;
  message << "Not implemented for solid: "
          << GetEntityType() << " !"
          << "\nReturning infinite boundinx box.";
  G4Exception("G4VSolid::BoundingLimits()", "GeomMgt1001",
              JustWarning, message);

  pMin.set(-kInfinity,-kInfinity,-kInfinity);
  pMax.set( kInfinity, kInfinity, kInfinity);
}

//////////////////////////////////////////////////////////////////////////
//
// Get G4VisExtent - bounding box for graphics

G4VisExtent G4VSolid::GetExtent () const
{
  G4VisExtent extent;
  G4VoxelLimits voxelLimits;  // Defaults to "infinite" limits.
  G4AffineTransform affineTransform;
  G4double vmin, vmax;
  CalculateExtent(kXAxis,voxelLimits,affineTransform,vmin,vmax);
  extent.SetXmin (vmin);
  extent.SetXmax (vmax);
  CalculateExtent(kYAxis,voxelLimits,affineTransform,vmin,vmax);
  extent.SetYmin (vmin);
  extent.SetYmax (vmax);
  CalculateExtent(kZAxis,voxelLimits,affineTransform,vmin,vmax);
  extent.SetZmin (vmin);
  extent.SetZmax (vmax);
  return extent;
}

G4Polyhedron* G4VSolid::CreatePolyhedron () const
{
  return nullptr;
}

G4Polyhedron* G4VSolid::GetPolyhedron () const
{
  return nullptr;
}
