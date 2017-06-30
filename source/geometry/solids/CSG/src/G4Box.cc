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
// $Id: G4Box.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
//
// Implementation for G4Box class
//
//  24.06.98 - V.Grichine: insideEdge in DistanceToIn(p,v)
//  20.09.98 - V.Grichine: new algorithm of DistanceToIn(p,v)
//  07.05.00 - V.Grichine: d= DistanceToIn(p,v), if d<e/2, d=0
//  09.06.00 - V.Grichine: safety in DistanceToIn(p) against Inside(p)=kOutside
//             and information before exception in DistanceToOut(p,v,...)
//  15.11.00 - D.Williams, V.Grichine: bug fixed in CalculateExtent - change
//                                     algorithm for rotated vertices
//  23.08.16 - E.Tcherniaev: use G4BoundingEnvelope for CalculateExtent()
//  18.04.17 - E.Tcherniaev: complete revision, speed-up
// --------------------------------------------------------------------

#include "G4Box.hh"

#if !defined(G4GEOM_USE_UBOX)

#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths

G4Box::G4Box(const G4String& pName,
                   G4double pX,
                   G4double pY,
                   G4double pZ)
  : G4CSGSolid(pName), fDx(pX), fDy(pY), fDz(pZ)
{
  delta = 0.5*kCarTolerance;
  if (pX < 2*kCarTolerance ||
      pY < 2*kCarTolerance ||
      pZ < 2*kCarTolerance)  // limit to thickness of surfaces
  {
    std::ostringstream message;
    message << "Dimensions too small for Solid: " << GetName() << "!" << G4endl
            << "     hX, hY, hZ = " << pX << ", " << pY << ", " << pZ;
    G4Exception("G4Box::G4Box()", "GeomSolids0002", FatalException, message);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4Box::G4Box( __void__& a )
  : G4CSGSolid(a), fDx(0.), fDy(0.), fDz(0.), delta(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Box::~G4Box()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Box::G4Box(const G4Box& rhs)
  : G4CSGSolid(rhs), fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz), delta(rhs.delta)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Box& G4Box::operator = (const G4Box& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fDx = rhs.fDx;
   fDy = rhs.fDy;
   fDz = rhs.fDz;
   delta = rhs.delta;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
//  Set X dimension

void G4Box::SetXHalfLength(G4double dx)
{
  if(dx > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    fDx = dx;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension X too small for solid: " << GetName() << "!"
            << G4endl
            << "       hX = " << dx;
    G4Exception("G4Box::SetXHalfLength()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fRebuildPolyhedron = true;
} 

//////////////////////////////////////////////////////////////////////////
//
//  Set Y dimension

void G4Box::SetYHalfLength(G4double dy) 
{
  if(dy > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    fDy = dy;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension Y too small for solid: " << GetName() << "!\n"
            << "       hY = " << dy;
    G4Exception("G4Box::SetYHalfLength()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fRebuildPolyhedron = true;
} 

//////////////////////////////////////////////////////////////////////////
//
//  Set Z dimension

void G4Box::SetZHalfLength(G4double dz)
{
  if(dz > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    fDz = dz;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension Z too small for solid: " << GetName() << "!\n"
            << "       hZ = " << dz;
    G4Exception("G4Box::SetZHalfLength()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fRebuildPolyhedron = true;
} 

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Box::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Box::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin.set(-fDx,-fDy,-fDz);
  pMax.set( fDx, fDy, fDz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Box::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Box::CalculateExtent(const EAxis pAxis,
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

//////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Box::Inside(const G4ThreeVector& p) const
{
  G4double dist = std::max(std::max(
                  std::abs(p.x())-fDx,
                  std::abs(p.y())-fDy),
                  std::abs(p.z())-fDz);
  if (dist > delta) return kOutside;
  return (dist > -delta) ? kSurface : kInside;
}

//////////////////////////////////////////////////////////////////////////
//
// Detect the side(s) and return corresponding normal

G4ThreeVector G4Box::SurfaceNormal( const G4ThreeVector& p) const
{
  G4ThreeVector norm(0,0,0);
  G4double px = p.x();
  if (std::abs(std::abs(px) - fDx) <= delta) norm.setX(px < 0 ? -1. : 1.);
  G4double py = p.y();
  if (std::abs(std::abs(py) - fDy) <= delta) norm.setY(py < 0 ? -1. : 1.);
  G4double pz = p.z();
  if (std::abs(std::abs(pz) - fDz) <= delta) norm.setZ(pz < 0 ? -1. : 1.);

  G4double nside = norm.mag2(); // number of sides = magnitude squared
  if (nside == 1)
    return norm;
  else if (nside > 1)
    return norm.unit(); // edge or corner
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
    G4cout.precision(oldprc);
    G4Exception("G4Box::SurfaceNormal(p)", "GeomSolids1002",
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

G4ThreeVector G4Box::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  G4double distx = std::abs(p.x()) - fDx;
  G4double disty = std::abs(p.y()) - fDy;
  G4double distz = std::abs(p.z()) - fDz;

  if (distx >= disty && distx >= distz)
    return G4ThreeVector(std::copysign(1.,p.x()), 0., 0.);
  if (disty >= distx && disty >= distz)
    return G4ThreeVector(0., std::copysign(1.,p.y()), 0.);
  else
    return G4ThreeVector(0., 0., std::copysign(1.,p.z()));
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection
//

G4double G4Box::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v) const
{
  // Check if point is on the surface and traveling away
  //
  if ((std::abs(p.x()) - fDx) >= -delta && p.x()*v.x() >= 0) return kInfinity;
  if ((std::abs(p.y()) - fDy) >= -delta && p.y()*v.y() >= 0) return kInfinity;
  if ((std::abs(p.z()) - fDz) >= -delta && p.z()*v.z() >= 0) return kInfinity;

  // Find intersection
  //
  G4double invx = (v.x() == 0) ? DBL_MAX : -1./v.x();
  G4double dx = std::copysign(fDx,invx);
  G4double txmin = (p.x() - dx)*invx;
  G4double txmax = (p.x() + dx)*invx;

  G4double invy = (v.y() == 0) ? DBL_MAX : -1./v.y();
  G4double dy = std::copysign(fDy,invy);
  G4double tymin = std::max(txmin,(p.y() - dy)*invy);
  G4double tymax = std::min(txmax,(p.y() + dy)*invy);

  G4double invz = (v.z() == 0) ? DBL_MAX : -1./v.z();
  G4double dz = std::copysign(fDz,invz);
  G4double tmin = std::max(tymin,(p.z() - dz)*invz);
  G4double tmax = std::min(tymax,(p.z() + dz)*invz);

  if (tmax <= tmin + delta) return kInfinity; // touch or no hit
  return (tmin < delta) ? 0. : tmin;
}

//////////////////////////////////////////////////////////////////////////
// 
// Appoximate distance to box.
// Returns largest perpendicular distance to the closest x/y/z sides of
// the box, which is the most fast estimation of the shortest distance to box
// - If inside return 0

G4double G4Box::DistanceToIn(const G4ThreeVector& p) const
{
  G4double dist = std::max(std::max(
                  std::abs(p.x())-fDx,
                  std::abs(p.y())-fDy),
                  std::abs(p.z())-fDz);
  return (dist > 0) ? dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of the box from inside and
// find normal at exit point, if required
// - when leaving the surface, return 0

G4double G4Box::DistanceToOut( const G4ThreeVector& p,
                               const G4ThreeVector& v,
                               const G4bool calcNorm,
                               G4bool *validNorm, G4ThreeVector *n) const
{
  // Check if point is on the surface and traveling away
  //
  if ((std::abs(p.x()) - fDx) >= -delta && p.x()*v.x() > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set((p.x() < 0) ? -1. : 1., 0., 0.);
    }
    return 0.;
  }
  if ((std::abs(p.y()) - fDy) >= -delta && p.y()*v.y() > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0., (p.y() < 0) ? -1. : 1., 0.);
    }
    return 0.;
  }
  if ((std::abs(p.z()) - fDz) >= -delta && p.z()*v.z() > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0., 0., (p.z() < 0) ? -1. : 1.);
    }
    return 0.;
  }

  // Find intersection
  //
  G4double vx = v.x();
  G4double tx = (vx == 0) ? DBL_MAX : (std::copysign(fDx,vx) - p.x())/vx;

  G4double vy = v.y();
  G4double ty = (vy == 0) ? tx : (std::copysign(fDy,vy) - p.y())/vy;
  G4double txy = std::min(tx,ty);

  G4double vz = v.z();
  G4double tz = (vz == 0) ? txy : (std::copysign(fDz,vz) - p.z())/vz;
  G4double tmax = std::min(txy,tz);

  // Set normal, if required, and return distance
  //
  if (calcNorm)
  {
    *validNorm = true;
    G4double xnew = p.x() + tmax*vx;
    if (std::abs(std::abs(xnew) - fDx) <= delta)
    {
      n->set(std::copysign(1.,xnew), 0., 0.);
      return tmax;
    }
    G4double ynew = p.y() + tmax*vy;
    if (std::abs(std::abs(ynew) - fDy) <= delta)
    {
      n->set(0., std::copysign(1.,ynew), 0.);
      return tmax;
    }
    else
    {
      G4double znew = p.z() + tmax*vz;
      n->set(0., 0., std::copysign(1.,znew));
    }
  }
  return tmax;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - if outside return 0

G4double G4Box::DistanceToOut(const G4ThreeVector& p) const
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
    G4cout.precision(oldprc);
    G4Exception("G4Box::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
  }
#endif
  G4double dist = std::max(std::max(
                  std::abs(p.x())-fDx,
                  std::abs(p.y())-fDy),
                  std::abs(p.z())-fDz);
  return (dist < 0) ? -dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Box::GetEntityType() const
{
  return G4String("G4Box");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Box::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << "Solid type: G4Box\n"
     << "Parameters: \n"
     << "   half length X: " << fDx/mm << " mm \n"
     << "   half length Y: " << fDy/mm << " mm \n"
     << "   half length Z: " << fDz/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// Return a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface

G4ThreeVector G4Box::GetPointOnSurface() const
{
  G4double sxy = fDx*fDy, sxz = fDx*fDz, syz = fDy*fDz;
  G4double select = (sxy + sxz + syz)*G4UniformRand();

  if (select < sxy)
    return G4ThreeVector((2.*G4UniformRand() - 1.)*fDx,
                         (2.*G4UniformRand() - 1.)*fDy,
                         (select < 0.5*sxy) ? -fDz : fDz);

  if (select < sxy + sxz)
    return G4ThreeVector((2.*G4UniformRand() - 1.)*fDx,
                         (select < sxy + 0.5*sxz) ? -fDy : fDy,
                         (2.*G4UniformRand() - 1.)*fDz);
  else
    return G4ThreeVector((select < sxy + sxz + 0.5*syz) ? -fDx : fDx,
                         (2.*G4UniformRand() - 1.)*fDy,
                         (2.*G4UniformRand() - 1.)*fDz);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Box::Clone() const
{
  return new G4Box(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Box::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

G4VisExtent G4Box::GetExtent() const
{
  return G4VisExtent (-fDx, fDx, -fDy, fDy, -fDz, fDz);
}

G4Polyhedron* G4Box::CreatePolyhedron () const
{
  return new G4PolyhedronBox (fDx, fDy, fDz);
}
#endif
