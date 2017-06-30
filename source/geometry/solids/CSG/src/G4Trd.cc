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
// $Id: G4Trd.cc 104894 2017-06-26 13:30:00Z gcosmo $
//
//
// Implementation for G4Trd class
//
// History:
//
// 25.05.17 E.Tcherniaev: complete revision, speed-up
// 23.09.16 E.Tcherniaev: use G4BoundingEnvelope for CalculateExtent(),
//                      removed CreateRotatedVertices()
// 28.04.05 V.Grichine: new SurfaceNormal according to J. Apostolakis proposal 
// 26.04.05, V.Grichine, new SurfaceNoramal is default
// 07.12.04, V.Grichine, SurfaceNoramal with edges/vertices.
// 07.05.00, V.Grichine, in d = DistanceToIn(p,v), if d<0.5*kCarTolerance, d=0
//    ~1996, V.Grichine, 1st implementation based on old code of P.Kent
//
//////////////////////////////////////////////////////////////////////////////

#include "G4Trd.hh"

#if !defined(G4GEOM_USE_UTRD)

#include "G4GeomTools.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
// Constructor - set & check half widths

G4Trd::G4Trd(const G4String& pName,
                   G4double pdx1, G4double pdx2,
                   G4double pdy1, G4double pdy2,
                   G4double pdz)
  : G4CSGSolid(pName), halfCarTolerance(0.5*kCarTolerance),
    fDx1(pdx1), fDx2(pdx2), fDy1(pdy1), fDy2(pdy2), fDz(pdz)
{
  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency
//
G4Trd::G4Trd( __void__& a )
  : G4CSGSolid(a), halfCarTolerance(0.5*kCarTolerance),
    fDx1(1.), fDx2(1.), fDy1(1.), fDy2(1.), fDz(1.)
{
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Trd::~G4Trd()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Trd::G4Trd(const G4Trd& rhs)
  : G4CSGSolid(rhs), halfCarTolerance(rhs.halfCarTolerance),
    fDx1(rhs.fDx1), fDx2(rhs.fDx2),
    fDy1(rhs.fDy1), fDy2(rhs.fDy2), fDz(rhs.fDz)
{
  for (G4int i=0; i<4; ++i) { fPlanes[i] = rhs.fPlanes[i]; }
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Trd& G4Trd::operator = (const G4Trd& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   halfCarTolerance = rhs.halfCarTolerance;
   fDx1 = rhs.fDx1; fDx2 = rhs.fDx2;
   fDy1 = rhs.fDy1; fDy2 = rhs.fDy2;
   fDz = rhs.fDz;
   for (G4int i=0; i<4; ++i) { fPlanes[i] = rhs.fPlanes[i]; }

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - set and check half-widths

void G4Trd::SetAllParameters(G4double pdx1, G4double pdx2,
                             G4double pdy1, G4double pdy2, G4double pdz)
{
  // Reset data of the base class
  fCubicVolume = 0;
  fSurfaceArea = 0;
  fRebuildPolyhedron = true;

  // Set parameters
  fDx1 = pdx1; fDx2 = pdx2;
  fDy1 = pdy1; fDy2 = pdy2;
  fDz  = pdz;

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Check dimensions

void G4Trd::CheckParameters()
{
  G4double dmin = 2*kCarTolerance;
  if ((fDx1 < 0 || fDx2 < 0 || fDy1 < 0 || fDy2 < 0 || fDz < dmin) ||
      (fDx1 < dmin && fDx2 < dmin) ||
      (fDy1 < dmin && fDy2 < dmin))
  {
    std::ostringstream message;
    message << "Invalid (too small or negative) dimensions for Solid: "
            << GetName()
            << "\n  X - " << fDx1 << ", " << fDx2
            << "\n  Y - " << fDy1 << ", " << fDy2
            << "\n  Z - " << fDz;
    G4Exception("G4Trd::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Set side planes

void G4Trd::MakePlanes()
{
  G4double dx = fDx1 - fDx2;
  G4double dy = fDy1 - fDy2;
  G4double dz = 2*fDz;
  G4double magx = std::sqrt(dx*dx + dz*dz);
  G4double magy = std::sqrt(dy*dy + dz*dz);

  // Set -Y & +Y planes
  //
  fPlanes[0].a =  0.;
  fPlanes[0].b = -dz/magy;
  fPlanes[0].c =  dy/magy;
  fPlanes[0].d = fPlanes[0].b*fDy1 + fPlanes[0].c*fDz;

  fPlanes[1].a =  fPlanes[0].a;
  fPlanes[1].b = -fPlanes[0].b;
  fPlanes[1].c =  fPlanes[0].c;
  fPlanes[1].d =  fPlanes[0].d;

  // Set -X & +X planes
  //
  fPlanes[2].a = -dz/magx;
  fPlanes[2].b =  0.;
  fPlanes[2].c =  dx/magx;
  fPlanes[2].d = fPlanes[2].a*fDx1 + fPlanes[2].c*fDz;

  fPlanes[3].a = -fPlanes[2].a;
  fPlanes[3].b =  fPlanes[2].b;
  fPlanes[3].c =  fPlanes[2].c;
  fPlanes[3].d =  fPlanes[2].d;
}

//////////////////////////////////////////////////////////////////////////
//
// Get volume

G4double G4Trd::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    fCubicVolume = 2*fDz*( (fDx1+fDx2)*(fDy1+fDy2) +
                           (fDx2-fDx1)*(fDy2-fDy1)/3 );
  }
  return fCubicVolume;
}

//////////////////////////////////////////////////////////////////////////
//
// Get surface area

G4double G4Trd::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    fSurfaceArea =
      4*(fDx1*fDy1+fDx2*fDy2) +
      2*(fDy1+fDy2)*std::hypot(fDx1-fDx2,2*fDz) +
      2*(fDx1+fDx2)*std::hypot(fDy1-fDy2,2*fDz);
  }
  return fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification

void G4Trd::ComputeDimensions(       G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Trd::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double dx1 = GetXHalfLength1();
  G4double dx2 = GetXHalfLength2();
  G4double dy1 = GetYHalfLength1();
  G4double dy2 = GetYHalfLength2();
  G4double dz  = GetZHalfLength();

  G4double xmax = std::max(dx1,dx2);
  G4double ymax = std::max(dy1,dy2);
  pMin.set(-xmax,-ymax,-dz);
  pMax.set( xmax, ymax, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Trd::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Trd::CalculateExtent( const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                     G4double& pMin, G4double& pMax ) const
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

  // Set bounding envelope (benv) and calculate extent
  //
  G4double dx1 = GetXHalfLength1();
  G4double dx2 = GetXHalfLength2();
  G4double dy1 = GetYHalfLength1();
  G4double dy2 = GetYHalfLength2();
  G4double dz  = GetZHalfLength();

  G4ThreeVectorList baseA(4), baseB(4);
  baseA[0].set(-dx1,-dy1,-dz);
  baseA[1].set( dx1,-dy1,-dz);
  baseA[2].set( dx1, dy1,-dz);
  baseA[3].set(-dx1, dy1,-dz);
  baseB[0].set(-dx2,-dy2, dz);
  baseB[1].set( dx2,-dy2, dz);
  baseB[2].set( dx2, dy2, dz);
  baseB[3].set(-dx2, dy2, dz);

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Trd::Inside( const G4ThreeVector& p ) const
{  
  G4double dx = fPlanes[3].a*std::abs(p.x())+fPlanes[3].c*p.z()+fPlanes[3].d;
  G4double dy = fPlanes[1].b*std::abs(p.y())+fPlanes[1].c*p.z()+fPlanes[1].d;
  G4double dxy = std::max(dx,dy);

  G4double dz = std::abs(p.z())-fDz;
  G4double dist = std::max(dz,dxy);

  if (dist > halfCarTolerance) return kOutside;
  return (dist > -halfCarTolerance) ? kSurface : kInside;
}

//////////////////////////////////////////////////////////////////////////
//
// Determine side where point is, and return corresponding normal

G4ThreeVector G4Trd::SurfaceNormal( const G4ThreeVector& p ) const
{
  // Check Z faces
  //
  G4double nz = 0;
  G4double dz = std::abs(p.z()) - fDz;
  if (std::abs(dz) <= halfCarTolerance) nz = (p.z() < 0) ? -1 : 1;

  // Check Y faces
  //
  G4double ny = 0;
  G4double dy1 = fPlanes[0].b*p.y();
  G4double dy2 = fPlanes[0].c*p.z() + fPlanes[0].d;
  if (std::abs(dy2 + dy1) <= halfCarTolerance)
  {
    ny += fPlanes[0].b;
    nz += fPlanes[0].c;
  }
  if (std::abs(dy2 - dy1) <= halfCarTolerance)
  {
    ny += fPlanes[1].b;
    nz += fPlanes[1].c;
  }

  // Check X faces
  //
  G4double nx = 0;
  G4double dx1 = fPlanes[2].a*p.x();
  G4double dx2 = fPlanes[2].c*p.z() + fPlanes[2].d;
  if (std::abs(dx2 + dx1) <= halfCarTolerance)
  {
    nx += fPlanes[2].a;
    nz += fPlanes[2].c;
  }
  if (std::abs(dx2 - dx1) <= halfCarTolerance)
  {
    nx += fPlanes[3].a;
    nz += fPlanes[3].c;
  }

  // Return normal
  //
  G4int nsurf = nx*nx + ny*ny + nz*nz + 0.5;                  // get magnitude
  if (nsurf == 1)      return G4ThreeVector(nx,ny,nz);
  else if (nsurf != 0) return G4ThreeVector(nx,ny,nz).unit(); // edge or corner
  else
  {
    // Point is not on the surface
    //
#ifdef G4CSGDEBUG
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is not on surface (!?) of solid: "
            << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc) ;
    G4Exception("G4Trd::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Trd::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  G4double dist = -DBL_MAX;
  G4int iside = 0;
  for (G4int i=0; i<4; ++i)
  {
    G4double d = fPlanes[i].a*p.x() +
                 fPlanes[i].b*p.y() +
                 fPlanes[i].c*p.z() + fPlanes[i].d;
    if (d > dist) { dist = d; iside = i; }
  }

  G4double distz = std::abs(p.z()) - fDz;
  if (dist > distz)
    return G4ThreeVector(fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c);
  else
    return G4ThreeVector(0, 0, (p.z() < 0) ? -1 : 1);
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside
//  - return kInfinity if no intersection

G4double G4Trd::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v ) const
{
  // Z intersections
  //
  if ((std::abs(p.z()) - fDz) >= -halfCarTolerance && p.z()*v.z() >= 0)
    return kInfinity;
  G4double invz = (-v.z() == 0) ? DBL_MAX : -1./v.z();
  G4double dz = (invz < 0) ? fDz : -fDz; 
  G4double tzmin = (p.z() + dz)*invz;
  G4double tzmax = (p.z() - dz)*invz;

  // Y intersections
  //
  G4double tmin0 = tzmin, tmax0 = tzmax;
  G4double ya = fPlanes[0].b*v.y(), yb = fPlanes[0].c*v.z();
  G4double yc = fPlanes[0].b*p.y(), yd = fPlanes[0].c*p.z()+fPlanes[0].d;
  G4double cos0 = yb + ya;
  G4double dis0 = yd + yc;
  if (dis0 >= -halfCarTolerance)
  {
    if (cos0 >= 0) return kInfinity;
    G4double tmp  = -dis0/cos0;
    if (tmin0 < tmp) tmin0 = tmp;
  }
  else if (cos0 > 0)
  {
    G4double tmp  = -dis0/cos0;
    if (tmax0 > tmp) tmax0 = tmp;
  }

  G4double tmin1 = tmin0, tmax1 = tmax0;
  G4double cos1 = yb - ya;
  G4double dis1 = yd - yc;
  if (dis1 >= -halfCarTolerance)
  {
    if (cos1 >= 0) return kInfinity;
    G4double tmp  = -dis1/cos1;
    if (tmin1 < tmp) tmin1 = tmp;
  }
  else if (cos1 > 0)
  {
    G4double tmp  = -dis1/cos1;
    if (tmax1 > tmp) tmax1 = tmp;
  }

  // X intersections
  //
  G4double tmin2 = tmin1, tmax2 = tmax1;
  G4double xa = fPlanes[2].a*v.x(), xb = fPlanes[2].c*v.z();
  G4double xc = fPlanes[2].a*p.x(), xd = fPlanes[2].c*p.z()+fPlanes[2].d;
  G4double cos2 = xb + xa;
  G4double dis2 = xd + xc;
  if (dis2 >= -halfCarTolerance)
  {
    if (cos2 >= 0) return kInfinity;
    G4double tmp  = -dis2/cos2;
    if (tmin2 < tmp) tmin2 = tmp;
  }
  else if (cos2 > 0)
  {
    G4double tmp  = -dis2/cos2;
    if (tmax2 > tmp) tmax2 = tmp;
  }

  G4double tmin3 = tmin2, tmax3 = tmax2;
  G4double cos3 = xb - xa;
  G4double dis3 = xd - xc;
  if (dis3 >= -halfCarTolerance)
  {
    if (cos3 >= 0) return kInfinity;
    G4double tmp  = -dis3/cos3;
    if (tmin3 < tmp) tmin3 = tmp;
  }
  else if (cos3 > 0)
  {
    G4double tmp  = -dis3/cos3;
    if (tmax3 > tmp) tmax3 = tmp;
  }

  // Find distance
  //
  G4double tmin = tmin3, tmax = tmax3;
  if (tmax <= tmin + halfCarTolerance) return kInfinity; // touch or no hit
  return (tmin < halfCarTolerance ) ? 0. : tmin;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// This is the best fast estimation of the shortest distance to trap
// - returns 0 if point is inside

G4double G4Trd::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double dx = fPlanes[3].a*std::abs(p.x())+fPlanes[3].c*p.z()+fPlanes[3].d;
  G4double dy = fPlanes[1].b*std::abs(p.y())+fPlanes[1].c*p.z()+fPlanes[1].d;
  G4double dxy = std::max(dx,dy);

  G4double dz = std::abs(p.z())-fDz;
  G4double dist = std::max(dz,dxy);

  return (dist > 0) ? dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside and
// find normal at exit point, if required
// - when leaving the surface, return 0

G4double G4Trd::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                              const G4bool calcNorm,
                                    G4bool *validNorm, G4ThreeVector *n) const
{
  // Z intersections
  //
  if ((std::abs(p.z()) - fDz) >= -halfCarTolerance && p.z()*v.z() > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0, 0, (p.z() < 0) ? -1 : 1);
    }
    return 0;
  }
  G4double vz = v.z();
  G4double tmax = (vz == 0) ? DBL_MAX : (std::copysign(fDz,vz) - p.z())/vz;
  G4int iside = (vz < 0) ? -4 : -2; // little trick: (-4+3)=-1, (-2+3)=+1

  // Y intersections
  //
  G4int i = 0;
  for ( ; i<2; ++i)
  {
    G4double cosa = fPlanes[i].b*v.y() + fPlanes[i].c*v.z();
    if (cosa > 0)
    {
      G4double dist = fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d;
      if (dist >= -halfCarTolerance)
      {
        if (calcNorm)
        {
          *validNorm = true;
          n->set(0, fPlanes[i].b, fPlanes[i].c);
        }
        return 0;
      }
      G4double tmp = -dist/cosa;
      if (tmax > tmp) { tmax = tmp; iside = i; }
    }
  }

  // X intersections
  //
  for ( ; i<4; ++i)
  {
    G4double cosa = fPlanes[i].a*v.x()+fPlanes[i].c*v.z();
    if (cosa > 0)
    {
      G4double dist = fPlanes[i].a*p.x()+fPlanes[i].c*p.z()+fPlanes[i].d;
      if (dist >= -halfCarTolerance)
      {
        if (calcNorm)
        {
           *validNorm = true;
           n->set(fPlanes[i].a, fPlanes[i].b, fPlanes[i].c);
        }
        return 0;
      }
      G4double tmp = -dist/cosa;
      if (tmax > tmp) { tmax = tmp; iside = i; }
    }
  }

  // Set normal, if required, and return distance
  //
  if (calcNorm) 
  {
    *validNorm = true;
    if (iside < 0)
      n->set(0, 0, iside + 3); // (-4+3)=-1, (-2+3)=+1
    else
      n->set(fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c);
  }
  return tmax;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - returns 0 if point is outside

G4double G4Trd::DistanceToOut( const G4ThreeVector& p ) const
{
#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: " << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc) ;
    G4Exception("G4Trd::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
  }
#endif
  G4double dx = fPlanes[3].a*std::abs(p.x())+fPlanes[3].c*p.z()+fPlanes[3].d;
  G4double dy = fPlanes[1].b*std::abs(p.y())+fPlanes[1].c*p.z()+fPlanes[1].d;
  G4double dxy = std::max(dx,dy);

  G4double dz = std::abs(p.z())-fDz;
  G4double dist = std::max(dz,dxy);
  
  return (dist < 0) ? -dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Trd::GetEntityType() const
{
  return G4String("G4Trd");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Trd::Clone() const
{
  return new G4Trd(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Trd::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Trd\n"
     << " Parameters: \n"
     << "    half length X, surface -dZ: " << fDx1/mm << " mm \n"
     << "    half length X, surface +dZ: " << fDx2/mm << " mm \n"
     << "    half length Y, surface -dZ: " << fDy1/mm << " mm \n"
     << "    half length Y, surface +dZ: " << fDy2/mm << " mm \n"
     << "    half length Z             : " <<  fDz/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Return a point randomly and uniformly selected on the solid surface

G4ThreeVector G4Trd::GetPointOnSurface() const
{
  // Set vertices
  //
  G4ThreeVector pt[8];
  pt[0].set(-fDx1,-fDy1,-fDz);
  pt[1].set( fDx1,-fDy1,-fDz);
  pt[2].set(-fDx1, fDy1,-fDz);
  pt[3].set( fDx1, fDy1,-fDz);
  pt[4].set(-fDx2,-fDy2, fDz);
  pt[5].set( fDx2,-fDy2, fDz);
  pt[6].set(-fDx2, fDy2, fDz);
  pt[7].set( fDx2, fDy2, fDz);

  // Set faces  (-Z, -Y, +Y, -X, +X, +Z)
  //
  G4int iface [6][4] =
    { {0,1,3,2}, {0,4,5,1}, {2,3,7,6}, {0,2,6,4}, {1,5,7,3}, {4,6,7,5} };

  // Set areas
  // 
  G4double syz = (fDy1 + fDy2)*std::hypot(fDx1 - fDx2, 2*fDz);
  G4double sxz = (fDx1 + fDx2)*std::hypot(fDy1 - fDy2, 2*fDz);
  G4double sface[6] = { 4*fDx1*fDy1, syz, syz, sxz, sxz, 4*fDx2*fDy2 };
  for (G4int i=1; i<6; ++i) { sface[i] += sface[i-1]; }

  // Select face
  //
  G4double select = sface[5]*G4UniformRand();
  G4int k = 5;
  if (select <= sface[4]) k = 4;
  if (select <= sface[3]) k = 3;
  if (select <= sface[2]) k = 2;
  if (select <= sface[1]) k = 1;
  if (select <= sface[0]) k = 0;

  // Select sub-triangle
  //
  G4int i0 = iface[k][0];
  G4int i1 = iface[k][1];
  G4int i2 = iface[k][2];
  G4int i3 = iface[k][3];
  G4double s1 = G4GeomTools::TriangleAreaNormal(pt[i0],pt[i1],pt[i3]).mag();
  G4double s2 = G4GeomTools::TriangleAreaNormal(pt[i2],pt[i1],pt[i3]).mag();
  if ((s1+s2)*G4UniformRand() > s1) i0 = i2;

  // Generate point
  //
  G4double u = G4UniformRand();
  G4double v = G4UniformRand();
  if (u + v > 1.) { u = 1. - u; v = 1. - v; }
  return (1.-u-v)*pt[i0] + u*pt[i1] + v*pt[i3];
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Trd::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Trd::CreatePolyhedron () const
{
  return new G4PolyhedronTrd2 (fDx1, fDx2, fDy1, fDy2, fDz);
}

#endif
