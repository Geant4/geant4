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
// class G4Ellipsoid
//
// Implementation of G4Ellipsoid class
//
// 10.11.99 G.Horton-Smith: first writing, based on G4Sphere class
// 25.02.05 G.Guerrieri: Revised
// 15.12.19 E.Tcherniaev: Complete revision
// --------------------------------------------------------------------

#include "G4Ellipsoid.hh"

#if !(defined(G4GEOM_USE_UELLIPSOID) && defined(G4GEOM_USE_SYS_USOLIDS))

#include "globals.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"
#include "G4RandomTools.hh"
#include "G4QuickRand.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex  = G4MUTEX_INITIALIZER;
  G4Mutex lateralareaMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4Ellipsoid::G4Ellipsoid(const G4String& name,
                               G4double xSemiAxis,
                               G4double ySemiAxis,
                               G4double zSemiAxis,
                               G4double zBottomCut,
                               G4double zTopCut)
  : G4VSolid(name), fDx(xSemiAxis), fDy(ySemiAxis), fDz(zSemiAxis),
    fZBottomCut(zBottomCut), fZTopCut(zTopCut)
{
  CheckParameters();
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4Ellipsoid::G4Ellipsoid( __void__& a )
  : G4VSolid(a), fDx(0.), fDy(0.), fDz(0.), fZBottomCut(0.), fZTopCut(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Ellipsoid::~G4Ellipsoid()
{
  delete fpPolyhedron; fpPolyhedron = nullptr;
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Ellipsoid::G4Ellipsoid(const G4Ellipsoid& rhs)
  : G4VSolid(rhs),
   fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz),
   fZBottomCut(rhs.fZBottomCut), fZTopCut(rhs.fZTopCut),
   halfTolerance(rhs.halfTolerance),
   fXmax(rhs.fXmax), fYmax(rhs.fYmax),
   fRsph(rhs.fRsph), fR(rhs.fR),
   fSx(rhs.fSx), fSy(rhs.fSy), fSz(rhs.fSz),
   fZMidCut(rhs.fZMidCut), fZDimCut(rhs.fZDimCut),
   fQ1(rhs.fQ1), fQ2(rhs.fQ2),
   fCubicVolume(rhs.fCubicVolume),
   fSurfaceArea(rhs.fSurfaceArea),
   fLateralArea(rhs.fLateralArea)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Ellipsoid& G4Ellipsoid::operator = (const G4Ellipsoid& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   fDx = rhs.fDx;
   fDy = rhs.fDy;
   fDz = rhs.fDz;
   fZBottomCut = rhs.fZBottomCut;
   fZTopCut = rhs.fZTopCut;

   halfTolerance = rhs.halfTolerance;
   fXmax = rhs.fXmax;
   fYmax = rhs.fYmax;
   fRsph = rhs.fRsph;
   fR = rhs.fR;
   fSx = rhs.fSx;
   fSy = rhs.fSy;
   fSz = rhs.fSz;
   fZMidCut = rhs.fZMidCut;
   fZDimCut = rhs.fZDimCut;
   fQ1 = rhs.fQ1;
   fQ2 = rhs.fQ2;

   fCubicVolume = rhs.fCubicVolume;
   fSurfaceArea = rhs.fSurfaceArea;
   fLateralArea = rhs.fLateralArea;

   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = nullptr;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Check parameters and make precalculation

void G4Ellipsoid::CheckParameters()
{
  halfTolerance = 0.5 * kCarTolerance; // half tolerance
  G4double dmin = 2. * kCarTolerance;

  // Check dimensions
  //
  if (fDx < dmin || fDy < dmin || fDz < dmin)
  {
    std::ostringstream message;
    message << "Invalid (too small or negative) dimensions for Solid: "
            << GetName()  << "\n"
            << "  semi-axis x: " << fDx << "\n"
            << "  semi-axis y: " << fDy << "\n"
            << "  semi-axis z: " << fDz;
    G4Exception("G4Ellipsoid::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
  G4double A = fDx;
  G4double B = fDy;
  G4double C = fDz;

  // Check cuts
  //
  if (fZBottomCut == 0. && fZTopCut == 0.)
  {
    fZBottomCut = -C;
    fZTopCut = C;
  }
  if (fZBottomCut >= C || fZTopCut <= -C || fZBottomCut >= fZTopCut)
  {
    std::ostringstream message;
    message << "Invalid Z cuts for Solid: "
            << GetName() << "\n"
            << "  bottom cut: " << fZBottomCut << "\n"
            << "  top cut: " << fZTopCut;
    G4Exception("G4Ellipsoid::CheckParameters()", "GeomSolids0002",
                FatalException, message);

  }
  fZBottomCut = std::max(fZBottomCut, -C);
  fZTopCut = std::min(fZTopCut, C);

  // Set extent in x and y
  fXmax = A;
  fYmax = B;
  if (fZBottomCut > 0.)
  {
    G4double ratio = fZBottomCut / C;
    G4double scale = std::sqrt((1. - ratio) * (1 + ratio));
    fXmax *= scale;
    fYmax *= scale;
  }
  if (fZTopCut < 0.)
  {
    G4double ratio  = fZTopCut / C;
    G4double scale  = std::sqrt((1. - ratio) * (1 + ratio));
    fXmax *= scale;
    fYmax *= scale;
  }

  // Set scale factors
  fRsph = std::max(std::max(A, B), C); // bounding sphere
  fR    = std::min(std::min(A, B), C); // radius of sphere after scaling
  fSx   = fR / A; // X scale factor
  fSy   = fR / B; // Y scale factor
  fSz   = fR / C; // Z scale factor

  // Scaled cuts
  fZMidCut = 0.5 * (fZTopCut + fZBottomCut) * fSz; // middle position
  fZDimCut = 0.5 * (fZTopCut - fZBottomCut) * fSz; // half distance

  // Coefficients for approximation of distance: Q1 * (x^2 + y^2 + z^2) - Q2
  fQ1 = 0.5 / fR;
  fQ2 = 0.5 * fR + halfTolerance * halfTolerance * fQ1;

  fCubicVolume = 0.; // volume
  fSurfaceArea = 0.; // surface area
  fLateralArea = 0.; // lateral surface area
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification

void G4Ellipsoid::ComputeDimensions(G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Ellipsoid::BoundingLimits(G4ThreeVector& pMin,
                                 G4ThreeVector& pMax) const
{
  pMin.set(-fXmax,-fYmax, fZBottomCut);
  pMax.set( fXmax, fYmax, fZTopCut);
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limits

G4bool
G4Ellipsoid::CalculateExtent(const EAxis pAxis,
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
// Return position of point: inside/outside/on surface

EInside G4Ellipsoid::Inside(const G4ThreeVector& p) const
{
  G4double x     = p.x() * fSx;
  G4double y     = p.y() * fSy;
  G4double z     = p.z() * fSz;
  G4double rr    = x * x + y * y + z * z;
  G4double distZ = std::abs(z - fZMidCut) - fZDimCut;
  G4double distR = fQ1 * rr - fQ2;
  G4double dist  = std::max(distZ, distR);

  if (dist > halfTolerance) return kOutside;
  return (dist > -halfTolerance) ? kSurface : kInside;
}

//////////////////////////////////////////////////////////////////////////
//
// Return unit normal to surface at p

G4ThreeVector G4Ellipsoid::SurfaceNormal( const G4ThreeVector& p) const
{
  G4ThreeVector norm(0., 0., 0.);
  G4int nsurf = 0;

  // Check cuts
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double z = p.z() * fSz;
  G4double distZ = std::abs(z - fZMidCut) - fZDimCut;
  if (std::abs(distZ) <= halfTolerance)
  {
    norm.setZ(std::copysign(1., z - fZMidCut));
    ++nsurf;
  }

  // Check lateral surface
  G4double distR = fQ1*(x*x + y*y + z*z) - fQ2;
  if (std::abs(distR) <= halfTolerance)
  {
    // normal = (p.x/a^2, p.y/b^2, p.z/c^2)
    norm += G4ThreeVector(x*fSx, y*fSy, z*fSz).unit();
    ++nsurf;
  }

  // Return normal
  if (nsurf == 1) return norm;
  else if (nsurf > 1) return norm.unit(); // edge
  else
  {
#ifdef G4SPECSDEBUG
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is not on surface (!?) of solid: "
            << GetName() << "\n";
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc);
    G4Exception("G4Ellipsoid::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Find surface nearest to point and return corresponding normal.
// This method normally should not be called.

G4ThreeVector G4Ellipsoid::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  G4double x  = p.x() * fSx;
  G4double y  = p.y() * fSy;
  G4double z  = p.z() * fSz;
  G4double rr = x * x + y * y + z * z;
  G4double distZ = std::abs(z - fZMidCut) - fZDimCut;
  G4double distR = std::sqrt(rr) - fR;
  if  (distR > distZ && rr > 0.) // distR > distZ is correct!
    return G4ThreeVector(x*fSx, y*fSy, z*fSz).unit();
  else
    return G4ThreeVector(0., 0., std::copysign(1., z - fZMidCut));
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside along normalised vector

G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector& p,
                                   const G4ThreeVector& v) const
{
  G4double offset = 0.;
  G4ThreeVector pcur = p;

  // Check if point is flying away, relative to bounding box
  //
  G4double safex = std::abs(p.x()) - fXmax;
  G4double safey = std::abs(p.y()) - fYmax;
  G4double safet = p.z() - fZTopCut;
  G4double safeb = fZBottomCut - p.z();

  if (safex >= -halfTolerance && p.x() * v.x() >= 0.) return kInfinity;
  if (safey >= -halfTolerance && p.y() * v.y() >= 0.) return kInfinity;
  if (safet >= -halfTolerance && v.z() >= 0.) return kInfinity;
  if (safeb >= -halfTolerance && v.z() <= 0.) return kInfinity;

  // Relocate point, if required
  //
  G4double safe = std::max(std::max(std::max(safex, safey), safet), safeb);
  if (safe > 32. * fRsph)
  {
    offset = (1. - 1.e-08) * safe - 2. * fRsph;
    pcur += offset * v;
    G4double dist = DistanceToIn(pcur, v);
    return (dist == kInfinity) ? kInfinity : dist + offset;
  }

  // Scale ellipsoid to sphere
  //
  G4double px = pcur.x() * fSx;
  G4double py = pcur.y() * fSy;
  G4double pz = pcur.z() * fSz;
  G4double vx = v.x() * fSx;
  G4double vy = v.y() * fSy;
  G4double vz = v.z() * fSz;

  // Check if point is leaving the solid
  //
  G4double dzcut = fZDimCut;
  G4double pzcut = pz - fZMidCut;
  G4double distZ = std::abs(pzcut) - dzcut;
  if (distZ >= -halfTolerance && pzcut * vz >= 0.) return kInfinity;

  G4double rr = px * px + py * py + pz * pz;
  G4double pv = px * vx + py * vy + pz * vz;
  G4double distR = fQ1 * rr - fQ2;
  if (distR >= -halfTolerance && pv >= 0.) return kInfinity;

  G4double A = vx * vx + vy * vy + vz * vz;
  G4double B = pv;
  G4double C = rr - fR * fR;
  G4double D = B * B - A * C;
  // scratch^2 = R^2 - (R - halfTolerance)^2 = 2 * R * halfTolerance
  G4double EPS = A * A * fR * kCarTolerance; // discriminant at scratching
  if (D <= EPS) return kInfinity; // no intersection or scratching

  // Find intersection with Z planes
  //
  G4double invz  = (vz == 0) ? DBL_MAX : -1./vz;
  G4double dz    = std::copysign(dzcut, invz);
  G4double tzmin = (pzcut - dz) * invz;
  G4double tzmax = (pzcut + dz) * invz;

  // Find intersection with lateral surface
  //
  G4double tmp = -B - std::copysign(std::sqrt(D), B);
  G4double t1 = tmp / A;
  G4double t2 = C / tmp;
  G4double trmin = std::min(t1, t2);
  G4double trmax = std::max(t1, t2);

  // Return distance
  //
  G4double tmin = std::max(tzmin, trmin);
  G4double tmax = std::min(tzmax, trmax);

  if (tmax - tmin <= halfTolerance) return kInfinity; // touch or no hit
  return (tmin < halfTolerance) ? offset : tmin + offset;
}

//////////////////////////////////////////////////////////////////////////
//
// Estimate distance to surface from outside

G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector& p) const
{
  G4double px = p.x();
  G4double py = p.y();
  G4double pz = p.z();

  // Safety distance to bounding box
  G4double distX = std::abs(px) - fXmax;
  G4double distY = std::abs(py) - fYmax;
  G4double distZ = std::max(pz - fZTopCut, fZBottomCut - pz);
  G4double distB = std::max(std::max(distX, distY), distZ);

  // Safety distance to lateral surface
  G4double x = px * fSx;
  G4double y = py * fSy;
  G4double z = pz * fSz;
  G4double distR = std::sqrt(x*x + y*y + z*z) - fR;

  // Return safety to in
  G4double dist = std::max(distB, distR);
  return (dist < 0.) ? 0. : dist;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface from inside along normalised vector

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p,
                                    const G4ThreeVector& v,
                                    const G4bool calcNorm,
                                          G4bool* validNorm,
                                          G4ThreeVector* n  ) const
{
  // Check if point flying away relative to Z planes
  //
  G4double pz = p.z() * fSz;
  G4double vz = v.z() * fSz;
  G4double dzcut = fZDimCut;
  G4double pzcut = pz - fZMidCut;
  G4double distZ = std::abs(pzcut) - dzcut;
  if (distZ >= -halfTolerance && pzcut * vz > 0.)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0., 0., std::copysign(1., pzcut));
    }
    return 0.;
  }

  // Check if point is flying away relative to lateral surface
  //
  G4double px = p.x() * fSx;
  G4double py = p.y() * fSy;
  G4double vx = v.x() * fSx;
  G4double vy = v.y() * fSy;
  G4double rr = px * px + py * py + pz * pz;
  G4double pv = px * vx + py * vy + pz * vz;
  G4double distR = fQ1 * rr - fQ2;
  if (distR >= -halfTolerance && pv > 0.)
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n = G4ThreeVector(px*fSx, py*fSy, pz*fSz).unit();
    }
    return 0.;
  }

  // Just in case check if point is outside (normally it should never be)
  //
  if (std::max(distZ, distR) > halfTolerance)
  {
#ifdef G4SPECSDEBUG
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: "
            << GetName() << G4endl;
    message << "Position:  " << p << G4endl;;
    message << "Direction: " << v;
    G4cout.precision(oldprc);
    G4Exception("G4Ellipsoid::DistanceToOut(p,v)", "GeomSolids1002",
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
  G4double A  = vx * vx + vy * vy + vz * vz;
  G4double B  = pv;
  G4double C  = rr - fR * fR;
  G4double D  = B * B - A * C;
  // It is expected that the point is located inside the sphere, so
  // max term in the expression for discriminant is A * R^2 and
  // max calculation error can be derived as follows:
  // A * (1 + 2e) * R^2 * (1 + 2e) = A * R^2 + (4 * A * R^2 * e)
  G4double EPS = 4. * A * fR * fR * DBL_EPSILON; // calculation error

  if (D <= EPS) // no intersection
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n = G4ThreeVector(px*fSx, py*fSy, pz*fSz).unit();
    }
    return 0.;
  }

  // Find intersection with Z cuts
  //
  G4double tzmax = (vz == 0.) ? DBL_MAX : (std::copysign(dzcut, vz) - pzcut) / vz;

  // Find intersection with lateral surface
  //
  G4double tmp = -B - std::copysign(std::sqrt(D), B);
  G4double trmax = (tmp < 0.) ? C/tmp : tmp/A;

  // Find distance and set normal, if required
  //
  G4double tmax = std::min(tzmax, trmax);
  //if (tmax < halfTolerance) tmax = 0.;

  if (calcNorm)
  {
    *validNorm = true;
    if (tmax == tzmax)
    {
      G4double pznew = pz + tmax * vz;
      n->set(0., 0., (pznew > fZMidCut) ? 1. : -1.);
    }
    else
    {
      G4double nx = (px + tmax * vx) * fSx;
      G4double ny = (py + tmax * vy) * fSy;
      G4double nz = (pz + tmax * vz) * fSz;
      *n = G4ThreeVector(nx, ny, nz).unit();
    }
  }
  return tmax;
}

//////////////////////////////////////////////////////////////////////////
//
// Estimate distance to surface from inside

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p) const
{
  // Safety distance in z direction
  G4double distZ = std::min(fZTopCut - p.z(), p.z() - fZBottomCut);

  // Safety distance to lateral surface
  G4double x = p.x() * fSx;
  G4double y = p.y() * fSy;
  G4double z = p.z() * fSz;
  G4double distR = fR - std::sqrt(x*x + y*y + z*z);

  // Return safety to out
  G4double dist = std::min(distZ, distR);
  return (dist < 0.) ? 0. : dist;
}

//////////////////////////////////////////////////////////////////////////
//
// Return entity type

G4GeometryType G4Ellipsoid::GetEntityType() const
{
  return G4String("G4Ellipsoid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4Ellipsoid::Clone() const
{
  return new G4Ellipsoid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to output stream

std::ostream& G4Ellipsoid::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters: \n"
     << "    semi-axis x: " << GetDx()/mm << " mm \n"
     << "    semi-axis y: " << GetDy()/mm << " mm \n"
     << "    semi-axis z: " << GetDz()/mm << " mm \n"
     << "    lower cut in z: " << GetZBottomCut()/mm << " mm \n"
     << "    upper cut in z: " << GetZTopCut()/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Return volume

G4double G4Ellipsoid::GetCubicVolume()
{
  if (fCubicVolume == 0.)
  {
    G4double piAB_3 = CLHEP::pi * fDx * fDy / 3.;
    fCubicVolume = 4. * piAB_3 * fDz;
    if (fZBottomCut > -fDz)
    {
      G4double hbot = 1. + fZBottomCut / fDz;
      fCubicVolume -= piAB_3 * hbot * hbot * (2. * fDz - fZBottomCut);
    }
    if (fZTopCut < fDz)
    {
      G4double htop = 1. - fZTopCut / fDz;
      fCubicVolume -= piAB_3 * htop * htop * (2. * fDz + fZTopCut);
    }
  }
  return fCubicVolume;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate area of lateral surface

G4double G4Ellipsoid::LateralSurfaceArea() const
{
  constexpr G4int NPHI = 1000.;
  constexpr G4double dPhi = CLHEP::halfpi/NPHI;
  constexpr G4double eps = 4.*DBL_EPSILON;

  G4double aa = fDx*fDx;
  G4double bb = fDy*fDy;
  G4double cc = fDz*fDz;
  G4double ab = fDx*fDy;
  G4double cc_aa = cc/aa;
  G4double cc_bb = cc/bb;
  G4double zmax = std::min(fZTopCut, fDz);
  G4double zmin = std::max(fZBottomCut,-fDz);
  G4double zmax_c = zmax/fDz;
  G4double zmin_c = zmin/fDz;
  G4double area = 0.;

  if (aa == bb) // spheroid, use analytical expression
  {
    G4double k = fDz/fDx;
    G4double kk = k*k;
    if (kk < 1. - eps)
    {
      G4double invk = fDx/fDz;
      G4double root = std::sqrt(1. - kk);
      G4double tmax = zmax_c*root;
      G4double tmin = zmin_c*root;
      area = CLHEP::pi*ab*
        ((zmax_c*std::sqrt(kk + tmax*tmax) - zmin_c*std::sqrt(kk + tmin*tmin)) +
         (std::asinh(tmax*invk) - std::asinh(tmin*invk))*kk/root);
    }
    else if (kk > 1. + eps)
    {
      G4double invk = fDx/fDz;
      G4double root = std::sqrt(kk - 1.);
      G4double tmax = zmax_c*root;
      G4double tmin = zmin_c*root;
      area = CLHEP::pi*ab*
        ((zmax_c*std::sqrt(kk - tmax*tmax) - zmin_c*std::sqrt(kk - tmin*tmin)) +
         (std::asin(tmax*invk) - std::asin(tmin*invk))*kk/root);
    }
    else
    {
      area = CLHEP::twopi*fDx*(zmax - zmin);
    }
    return area;
  }

  // ellipsoid, integration along phi
  for (G4int i = 0; i < NPHI; ++i)
  {
    G4double sinPhi = std::sin(dPhi*(i + 0.5));
    G4double kk = cc_aa + (cc_bb - cc_aa)*sinPhi*sinPhi;
    if (kk < 1. - eps)
    {
      G4double root = std::sqrt(1. - kk);
      G4double tmax = zmax_c*root;
      G4double tmin = zmin_c*root;
      G4double invk = 1./std::sqrt(kk);
      area += 2.*ab*dPhi*
        ((zmax_c*std::sqrt(kk + tmax*tmax) - zmin_c*std::sqrt(kk + tmin*tmin)) +
         (std::asinh(tmax*invk) - std::asinh(tmin*invk))*kk/root);
    }
    else if (kk > 1. + eps)
    {
      G4double root = std::sqrt(kk - 1.);
      G4double tmax = zmax_c*root;
      G4double tmin = zmin_c*root;
      G4double invk = 1./std::sqrt(kk);
      area += 2.*ab*dPhi*
        ((zmax_c*std::sqrt(kk - tmax*tmax) - zmin_c*std::sqrt(kk - tmin*tmin)) +
         (std::asin(tmax*invk) - std::asin(tmin*invk))*kk/root);
    }
    else
    {
      area += 4.*ab*dPhi*(zmax_c - zmin_c);
    }
  }
  return area;
}

//////////////////////////////////////////////////////////////////////////
//
// Return surface area

G4double G4Ellipsoid::GetSurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    G4double piAB = CLHEP::pi * fDx * fDy;
    fSurfaceArea = LateralSurfaceArea();
    if (fZBottomCut > -fDz)
    {
      G4double hbot = 1. + fZBottomCut / fDz;
      fSurfaceArea += piAB * hbot * (2. - hbot);
    }
    if (fZTopCut < fDz)
    {
      G4double htop = 1. - fZTopCut / fDz;
      fSurfaceArea += piAB * htop * (2. - htop);
    }
  }
  return fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Return random point on surface

G4ThreeVector G4Ellipsoid::GetPointOnSurface() const
{
  G4double A    = GetDx();
  G4double B    = GetDy();
  G4double C    = GetDz();
  G4double Zbot = GetZBottomCut();
  G4double Ztop = GetZTopCut();

  // Calculate cut areas
  G4double Hbot = 1. + Zbot / C;
  G4double Htop = 1. - Ztop / C;
  G4double piAB = CLHEP::pi * A * B;
  G4double Sbot = piAB * Hbot * (2. - Hbot);
  G4double Stop = piAB * Htop * (2. - Htop);

  // Get area of lateral surface
  if (fLateralArea == 0.)
  {
    G4AutoLock l(&lateralareaMutex);
    fLateralArea = LateralSurfaceArea();
    l.unlock();
  }
  G4double Slat = fLateralArea;

  // Select surface (0 - bottom cut, 1 - lateral surface, 2 - top cut)
  G4double select = (Sbot + Slat + Stop) * G4QuickRand();
  G4int k = 0;
  if (select > Sbot) k = 1;
  if (select > Sbot + Slat) k = 2;

  // Pick random point on selected surface (rejection sampling)
  G4ThreeVector p;
  switch (k)
  {
    case 0: // bootom z-cut
    {
      G4double scale = std::sqrt(Hbot * (2. - Hbot));
      G4TwoVector rho = G4RandomPointInEllipse(A * scale, B * scale);
      p.set(rho.x(), rho.y(), Zbot);
      break;
    }
    case 1: // lateral surface
    {
      G4double x, y, z;
      G4double mu_max = std::max(std::max(A * B, A * C), B * C);
      for (G4int i = 0; i < 1000; ++i)
      {
        // generate random point on unit sphere
        z = (Zbot + (Ztop - Zbot) * G4QuickRand()) / C;
        G4double rho = std::sqrt((1. + z) * (1. - z));
        G4double phi = CLHEP::twopi * G4QuickRand();
        x = rho * std::cos(phi);
        y = rho * std::sin(phi);
        // check  acceptance
        G4double xbc = x * B * C;
        G4double yac = y * A * C;
        G4double zab = z * A * B;
        G4double mu  = std::sqrt(xbc * xbc + yac * yac + zab * zab);
        if (mu_max * G4QuickRand() <= mu) break;
      }
      p.set(A * x, B * y, C * z);
      break;
    }
    case 2: // top z-cut
    {
      G4double scale  = std::sqrt(Htop * (2. - Htop));
      G4TwoVector rho = G4RandomPointInEllipse(A * scale, B * scale);
      p.set(rho.x(), rho.y(), Ztop);
      break;
    }
  }
  return p;
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Ellipsoid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Return vis extent

G4VisExtent G4Ellipsoid::GetExtent() const
{
  return G4VisExtent(-fXmax, fXmax, -fYmax, fYmax, fZBottomCut, fZTopCut);
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron

G4Polyhedron* G4Ellipsoid::CreatePolyhedron () const
{
  return new G4PolyhedronEllipsoid(fDx, fDy, fDz, fZBottomCut, fZTopCut);
}

//////////////////////////////////////////////////////////////////////////
//
// Return pointer to polyhedron

G4Polyhedron* G4Ellipsoid::GetPolyhedron () const
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

#endif // !defined(G4GEOM_USE_UELLIPSOID) || !defined(G4GEOM_USE_SYS_USOLIDS)
