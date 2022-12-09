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
// G4EllipticalTube implementation
//
// Author: David C. Williams (davidw@scipp.ucsc.edu)
// Revision: Evgueni Tcherniaev (evgueni.tcherniaev@cern.ch), 23.12.2019
// --------------------------------------------------------------------

#include "G4EllipticalTube.hh"

#if !(defined(G4GEOM_USE_UELLIPTICALTUBE) && defined(G4GEOM_USE_SYS_USOLIDS))

#include "G4GeomTools.hh"
#include "G4RandomTools.hh"
#include "G4ClippablePolygon.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4BoundingEnvelope.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4EllipticalTube::G4EllipticalTube( const G4String &name,
                                          G4double Dx,
                                          G4double Dy,
                                          G4double Dz )
  : G4VSolid(name), fDx(Dx), fDy(Dy), fDz(Dz)
{
  CheckParameters();
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4EllipticalTube::G4EllipticalTube( __void__& a )
  : G4VSolid(a), halfTolerance(0.), fDx(0.), fDy(0.), fDz(0.),
    fRsph(0.), fDDx(0.), fDDy(0.), fSx(0.), fSy(0.), fR(0.),
    fQ1(0.), fQ2(0.), fScratch(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4EllipticalTube::~G4EllipticalTube()
{
  delete fpPolyhedron; fpPolyhedron = nullptr;
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4EllipticalTube::G4EllipticalTube(const G4EllipticalTube& rhs)
  : G4VSolid(rhs), halfTolerance(rhs.halfTolerance),
    fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    fRsph(rhs.fRsph), fDDx(rhs.fDDx), fDDy(rhs.fDDy),
    fSx(rhs.fSx), fSy(rhs.fSy), fR(rhs.fR),
    fQ1(rhs.fQ1), fQ2(rhs.fQ2), fScratch(rhs.fScratch)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4EllipticalTube& G4EllipticalTube::operator = (const G4EllipticalTube& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   halfTolerance = rhs.halfTolerance;
   fDx = rhs.fDx;
   fDy = rhs.fDy;
   fDz = rhs.fDz;
   fCubicVolume = rhs.fCubicVolume;
   fSurfaceArea = rhs.fSurfaceArea;

   fRsph = rhs.fRsph;
   fDDx  = rhs.fDDx;
   fDDy  = rhs.fDDy;
   fSx   = rhs.fSx;
   fSy   = rhs.fSy;
   fR    = rhs.fR;
   fQ1   = rhs.fQ1;
   fQ2   = rhs.fQ2;
   fScratch = rhs.fScratch;

   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = nullptr;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Check dimensions

void G4EllipticalTube::CheckParameters()
{
  // Check dimensions
  //
  halfTolerance = 0.5*kCarTolerance; // half tolerance
  G4double dmin = 2*kCarTolerance;
  if (fDx < dmin || fDy < dmin || fDz < dmin)
  {
    std::ostringstream message;
    message << "Invalid (too small or negative) dimensions for Solid: "
            << GetName()
            << "\n  Dx = " << fDx
            << "\n  Dy = " << fDy
            << "\n  Dz = " << fDz;
    G4Exception("G4EllipticalTube::CheckParameters()", "GeomSolids0002",
	        FatalException, message);
  }

  // Set pre-calculatated values
  //
  halfTolerance = 0.5*kCarTolerance; // half tolerance
  fRsph = std::sqrt(fDx * fDx + fDy * fDy + fDz * fDz); // radius of surrounding sphere
  fDDx = fDx * fDx; // X semi-axis squared
  fDDy = fDy * fDy; // Y semi-axis squared

  fR = std::min(fDx, fDy); // resulting radius, after scaling elipse to circle
  fSx = fR / fDx; // X scale factor
  fSy = fR / fDy; // Y scale factor

  fQ1 = 0.5 / fR; // distance approxiamtion dist = Q1 * (x^2 + y^2) - Q2
  fQ2 = 0.5 * (fR + halfTolerance * halfTolerance / fR);
  fScratch = 2. * fR * fR * DBL_EPSILON; // scratch within calculation error thickness
  // fScratch = (B * B / A) * (2. + halfTolerance / A) * halfTolerance; // alternative
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4EllipticalTube::BoundingLimits( G4ThreeVector& pMin,
                                       G4ThreeVector& pMax ) const
{
  pMin.set(-fDx,-fDy,-fDz);
  pMax.set( fDx, fDy, fDz);
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4EllipticalTube::CalculateExtent( const EAxis pAxis,
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
  return bbox.CalculateExtent(pAxis,pVoxelLimit, pTransform, pMin, pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis, pVoxelLimit, pTransform, pMin, pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  G4double dx = fDx;
  G4double dy = fDy;
  G4double dz = fDz;

  // Set bounding envelope (benv) and calculate extent
  //
  const G4int NSTEPS = 24; // number of steps for whole circle
  G4double ang = twopi/NSTEPS;

  G4double sinHalf = std::sin(0.5*ang);
  G4double cosHalf = std::cos(0.5*ang);
  G4double sinStep = 2.*sinHalf*cosHalf;
  G4double cosStep = 1. - 2.*sinHalf*sinHalf;
  G4double sx = dx/cosHalf;
  G4double sy = dy/cosHalf;

  G4double sinCur = sinHalf;
  G4double cosCur = cosHalf;
  G4ThreeVectorList baseA(NSTEPS),baseB(NSTEPS);
  for (G4int k=0; k<NSTEPS; ++k)
  {
    baseA[k].set(sx*cosCur,sy*sinCur,-dz);
    baseB[k].set(sx*cosCur,sy*sinCur, dz);

    G4double sinTmp = sinCur;
    sinCur = sinCur*cosStep + cosCur*sinStep;
    cosCur = cosCur*cosStep - sinTmp*sinStep;
  }

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;
  G4BoundingEnvelope benv(bmin, bmax, polygons);
  exist = benv.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Determine where is point: inside, outside or on surface
//

EInside G4EllipticalTube::Inside( const G4ThreeVector& p ) const
{
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double distR = fQ1 * (x * x + y * y) - fQ2;
  G4double distZ = std::abs(p.z()) - fDz;
  G4double dist = std::max(distR, distZ);

  if (dist > halfTolerance) return kOutside;
  return (dist > -halfTolerance) ? kSurface : kInside;
}

//////////////////////////////////////////////////////////////////////////
//
// Return unit normal at surface closest to p

G4ThreeVector G4EllipticalTube::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4ThreeVector norm(0, 0, 0);
  G4int nsurf = 0;

  // check lateral surface
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double distR = fQ1 * (x * x + y * y) - fQ2;
  if (std::abs(distR) <= halfTolerance)
  {
    norm = G4ThreeVector(p.x() * fDDy, p.y() * fDDx, 0.).unit();
    ++nsurf;
  }

  // check lateral bases
  G4double distZ = std::abs(p.z()) - fDz;
  if (std::abs(distZ) <= halfTolerance)
  {
    norm.setZ(p.z() < 0 ? -1. : 1.);
    ++nsurf;
  }

  // return normal
  if (nsurf == 1) return norm;
  else if (nsurf > 1) return norm.unit(); // edge
  else
  {
    // Point is not on the surface
    //
#ifdef G4SPECDEBUG
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is not on surface (!?) of solid: "
            << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc);
    G4Exception("G4EllipticalTube::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Find surface nearest to point and return corresponding normal.
// The algorithm is similar to the algorithm used in Inside().
// This method normally should not be called.

G4ThreeVector
G4EllipticalTube::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double distR = fQ1 * (x * x + y * y) - fQ2;
  G4double distZ = std::abs(p.z()) - fDz;
  if (distR > distZ && (x * x + y * y) > 0)
    return G4ThreeVector(p.x() * fDDy, p.y() * fDDx, 0.).unit();
  else
    return G4ThreeVector(0, 0, (p.z() < 0 ? -1. : 1.));
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector,
// return kInfinity if no intersection, or distance < halfTolerance

G4double G4EllipticalTube::DistanceToIn( const G4ThreeVector& p,
                                         const G4ThreeVector& v ) const
{
  G4double offset = 0.;
  G4ThreeVector pcur = p;

  // Check if point is flying away
  //
  G4double safex = std::abs(pcur.x()) - fDx;
  G4double safey = std::abs(pcur.y()) - fDy;
  G4double safez = std::abs(pcur.z()) - fDz;

  if (safez >= -halfTolerance && pcur.z() * v.z() >= 0.) return kInfinity;
  if (safey >= -halfTolerance && pcur.y() * v.y() >= 0.) return kInfinity;
  if (safex >= -halfTolerance && pcur.x() * v.x() >= 0.) return kInfinity;

  // Relocate point, if required
  //
  G4double Dmax = 32. * fRsph;
  if (std::max(std::max(safex, safey), safez) > Dmax)
  {
    offset = (1. - 1.e-08) * pcur.mag() - 2. * fRsph;
    pcur += offset * v;
    G4double dist = DistanceToIn(pcur, v);
    return (dist == kInfinity) ? kInfinity : dist + offset;
  }

  // Scale elliptical tube to cylinder
  //
  G4double px = pcur.x() * fSx;
  G4double py = pcur.y() * fSy;
  G4double pz = pcur.z();
  G4double vx = v.x() * fSx;
  G4double vy = v.y() * fSy;
  G4double vz = v.z();

  // Set coefficients of quadratic equation: A t^2 + 2B t + C = 0
  //
  G4double rr = px * px + py * py;
  G4double A  = vx * vx + vy * vy;
  G4double B  = px * vx + py * vy;
  G4double C  = rr - fR * fR;
  G4double D  = B * B - A * C;

  // Check if point is flying away relative to lateral surface
  //
  G4double distR  = fQ1 * rr - fQ2;
  G4bool parallelToZ = (A < DBL_EPSILON || std::abs(vz) >= 1.);
  if (distR >= -halfTolerance && (B >= 0. || parallelToZ)) return kInfinity;

  // Find intersection with Z planes
  //
  G4double invz  = (vz == 0) ? DBL_MAX : -1./vz;
  G4double dz    = std::copysign(fDz, invz);
  G4double tzmin = (pz - dz) * invz;
  G4double tzmax = (pz + dz) * invz;

  // Solve qudratic equation. There are two cases special where D <= 0:
  //   1) trajectory parallel to Z axis (A = 0, B = 0, C - any, D = 0)
  //   2) touch (D = 0) or no intersection (D < 0) with lateral surface
  //
  if (parallelToZ) return (tzmin<halfTolerance) ? offset : tzmin + offset; // 1)
  if (D <= A * A * fScratch) return kInfinity; // 2)

  // Find roots of quadratic equation
  G4double tmp = -B - std::copysign(std::sqrt(D), B);
  G4double t1 = tmp / A;
  G4double t2 = C / tmp;
  G4double trmin = std::min(t1, t2);
  G4double trmax = std::max(t1, t2);

  // Return distance
  G4double tin  = std::max(tzmin, trmin);
  G4double tout = std::min(tzmax, trmax);

  if (tout <= tin + halfTolerance) return kInfinity; // touch or no hit
  return (tin<halfTolerance) ? offset : tin + offset;
}

//////////////////////////////////////////////////////////////////////////
//
// Estimate distance to the surface from outside,
// returns 0 if point is inside

G4double G4EllipticalTube::DistanceToIn( const G4ThreeVector& p ) const
{
  // safety distance to bounding box
  G4double distX = std::abs(p.x()) - fDx;
  G4double distY = std::abs(p.y()) - fDy;
  G4double distZ = std::abs(p.z()) - fDz;
  G4double distB = std::max(std::max(distX, distY), distZ);
  // return (distB < 0) ? 0 : distB;

  // safety distance to lateral surface
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double distR = std::sqrt(x * x + y * y) - fR;

  // return SafetyToIn
  G4double dist = std::max(distB, distR);
  return (dist < 0) ? 0 : dist;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from inside and find normal
// at exit point, if required
// - when leaving the surface, return 0

G4double G4EllipticalTube::DistanceToOut( const G4ThreeVector& p,
                                          const G4ThreeVector& v,
                                          const G4bool calcNorm,
                                                G4bool* validNorm,
                                                G4ThreeVector* n ) const
{
  // Check if point flying away relative to Z planes
  //
  G4double pz = p.z();
  G4double vz = v.z();
  G4double distZ = std::abs(pz) - fDz;
  if (distZ >= -halfTolerance && pz * vz > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0, 0, (pz < 0) ? -1. : 1.);
    }
    return 0.;
  }
  G4double tzmax = (vz == 0) ? DBL_MAX : (std::copysign(fDz, vz) - pz) / vz;

  // Scale elliptical tube to cylinder
  //
  G4double px = p.x() * fSx;
  G4double py = p.y() * fSy;
  G4double vx = v.x() * fSx;
  G4double vy = v.y() * fSy;

  // Check if point is flying away relative to lateral surface
  //
  G4double rr = px * px + py * py;
  G4double B  = px * vx + py * vy;
  G4double distR  = fQ1 * rr - fQ2;
  if (distR >= -halfTolerance && B > 0.)
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n = G4ThreeVector(px * fDDy, py * fDDx, 0.).unit();
    }
    return 0.;
  }

  // Just in case check if point is outside, normally it should never be
  //
  if (std::max(distZ, distR) > halfTolerance)
  {
#ifdef G4SPECDEBUG
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: "
            << GetName() << G4endl;
    message << "Position:  " << p << G4endl;;
    message << "Direction: " << v;
    G4cout.precision(oldprc);
    G4Exception("G4EllipticalTube::DistanceToOut(p,v)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    if (calcNorm)
    {
      *validNorm = true;
      *n = ApproxSurfaceNormal(p);
    }
    return 0.;
  }

  // Set coefficients of quadratic equation: A t^2 + 2B t + C = 0
  //
  G4double A  = vx * vx + vy * vy;
  G4double C  = rr - fR * fR;
  G4double D  = B * B - A * C;

  // Solve qudratic equation. There are two special cases where D <= 0:
  //   1) trajectory parallel to Z axis (A = 0, B = 0, C - any, D = 0)
  //   2) touch (D = 0) or no intersection (D < 0) with lateral surface
  //
  G4bool parallelToZ = (A < DBL_EPSILON || std::abs(vz) >= 1.);
  if (parallelToZ) // 1)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0, 0, (vz < 0) ? -1. : 1.);
    }
    return tzmax;
  }
  if (D <= A * A * fScratch) // 2)
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n = G4ThreeVector(px * fDDy, py * fDDx, 0.).unit();
    }
    return 0.;
  }

  // Find roots of quadratic equation
  G4double tmp = -B - std::copysign(std::sqrt(D), B);
  G4double t1 = tmp / A;
  G4double t2 = C / tmp;
  G4double trmax = std::max(t1, t2);

  // Return distance
  G4double tmax = std::min(tzmax, trmax);

  // Set normal, if required, and return distance
  //
  if (calcNorm)
  {
    *validNorm = true;
    G4ThreeVector pnew = p + tmax * v;
    if (tmax == tzmax)
      n->set(0, 0, (pnew.z() < 0) ? -1. : 1.);
    else
      *n = G4ThreeVector(pnew.x() * fDDy, pnew.y() * fDDx, 0.).unit();
  }
  return tmax;
}

//////////////////////////////////////////////////////////////////////////
//
// Estimate distance to the surface from inside,
// returns 0 if point is outside
//

G4double G4EllipticalTube::DistanceToOut( const G4ThreeVector& p ) const
{
#ifdef G4SPECDEBUG
  if( Inside(p) == kOutside )
  {
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: " << GetName() << "\n"
            << "Position:\n"
            << "   p.x() = "  << p.x()/mm << " mm\n"
            << "   p.y() = "  << p.y()/mm << " mm\n"
            << "   p.z() = "  << p.z()/mm << " mm";
    message.precision(oldprc) ;
    G4Exception("G4ElliptocalTube::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, message);
    DumpInfo();
  }
#endif
  // safety distance to Z-bases
  G4double distZ = fDz - std::abs(p.z());

  // safety distance lateral surface
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double distR = fR - std::sqrt(x * x + y * y);

  // return SafetyToOut
  G4double dist = std::min(distZ, distR);
  return (dist < 0) ? 0 : dist;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4EllipticalTube::GetEntityType() const
{
  return G4String("G4EllipticalTube");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4EllipticalTube::Clone() const
{
  return new G4EllipticalTube(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Return volume

G4double G4EllipticalTube::GetCubicVolume()
{
  if (fCubicVolume == 0.)
  {
    fCubicVolume = twopi * fDx * fDy * fDz;
  }
  return fCubicVolume;
}

//////////////////////////////////////////////////////////////////////////
//
// Return cached surface area

G4double G4EllipticalTube::GetCachedSurfaceArea() const
{
  G4ThreadLocalStatic G4double cached_Dx = 0;
  G4ThreadLocalStatic G4double cached_Dy = 0;
  G4ThreadLocalStatic G4double cached_Dz = 0;
  G4ThreadLocalStatic G4double cached_area = 0;
  if (cached_Dx != fDx || cached_Dy != fDy || cached_Dz != fDz)
  {
    cached_Dx = fDx;
    cached_Dy = fDy;
    cached_Dz = fDz;
    cached_area = 2.*(pi*fDx*fDy + G4GeomTools::EllipsePerimeter(fDx, fDy)*fDz);
  }
  return cached_area;
}

//////////////////////////////////////////////////////////////////////////
//
// Return surface area

G4double G4EllipticalTube::GetSurfaceArea()
{
  if(fSurfaceArea == 0.)
  {
    fSurfaceArea = GetCachedSurfaceArea();
  }
  return fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to output stream

std::ostream& G4EllipticalTube::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4EllipticalTube\n"
     << " Parameters: \n"
     << "    length Z: " << fDz/mm << " mm \n"
     << "    lateral surface equation: \n"
     << "       (X / " << fDx << ")^2 + (Y / " << fDy << ")^2 = 1 \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


//////////////////////////////////////////////////////////////////////////
//
// Pick up a random point on the surface 

G4ThreeVector G4EllipticalTube::GetPointOnSurface() const
{
  // Select surface (0 - base at -Z, 1 - base at +Z, 2 - lateral surface)
  //
  G4double sbase = pi * fDx * fDy;
  G4double ssurf = GetCachedSurfaceArea();
  G4double select = ssurf * G4UniformRand();

  G4int k = 0;
  if (select > sbase) k = 1;
  if (select > 2. * sbase) k = 2;

  // Pick random point on selected surface (rejection sampling)
  //
  G4ThreeVector p;
  switch (k) {
    case 0: // base at -Z
    {
      G4TwoVector rho = G4RandomPointInEllipse(fDx, fDy);
      p.set(rho.x(), rho.y(), -fDz);
      break;
    }
    case 1: // base at +Z
    {
      G4TwoVector rho = G4RandomPointInEllipse(fDx, fDy);
      p.set(rho.x(), rho.y(), fDz);
      break;
    }
    case 2: // lateral surface
    {
      G4TwoVector rho = G4RandomPointOnEllipse(fDx, fDy);
      p.set(rho.x(), rho.y(), (2. * G4UniformRand() - 1.) * fDz);
      break;
    }
  }
  return p;
}


//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron

G4Polyhedron* G4EllipticalTube::CreatePolyhedron() const
{
  // create cylinder with radius=1...
  //
  G4Polyhedron* eTube = new G4PolyhedronTube(0., 1., fDz);

  // apply non-uniform scaling...
  //
  eTube->Transform(G4Scale3D(fDx, fDy, 1.));
  return eTube;
}

//////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron

G4Polyhedron* G4EllipticalTube::GetPolyhedron () const
{
  if (fpPolyhedron == nullptr ||
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

//////////////////////////////////////////////////////////////////////////
//
// DescribeYourselfTo

void G4EllipticalTube::DescribeYourselfTo( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

//////////////////////////////////////////////////////////////////////////
//
// GetExtent

G4VisExtent G4EllipticalTube::GetExtent() const
{
  return G4VisExtent( -fDx, fDx, -fDy, fDy, -fDz, fDz );
}

#endif // !defined(G4GEOM_USE_UELLIPTICALTUBE) || !defined(G4GEOM_USE_SYS_USOLIDS)
