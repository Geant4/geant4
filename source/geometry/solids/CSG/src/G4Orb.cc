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
// Implementation for G4Orb class
//
// 20.08.03 V.Grichine - created
// 08.08.17 E.Tcherniaev - complete revision, speed-up
// --------------------------------------------------------------------

#include "G4Orb.hh"

#if !defined(G4GEOM_USE_UORB)

#include "G4TwoVector.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "G4RandomDirection.hh"
#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4Orb::G4Orb( const G4String& pName, G4double pRmax )
  : G4CSGSolid(pName), fRmax(pRmax)
{
  Initialize();
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4Orb::G4Orb( __void__& a )
  : G4CSGSolid(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Orb::~G4Orb()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Orb::G4Orb(const G4Orb& rhs)
  : G4CSGSolid(rhs), fRmax(rhs.fRmax), halfRmaxTol(rhs.halfRmaxTol),
    sqrRmaxPlusTol(rhs.sqrRmaxPlusTol), sqrRmaxMinusTol(rhs.sqrRmaxMinusTol)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Orb& G4Orb::operator = (const G4Orb& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fRmax = rhs.fRmax;
   halfRmaxTol = rhs.halfRmaxTol;
   sqrRmaxPlusTol = rhs.sqrRmaxPlusTol;
   sqrRmaxMinusTol = rhs.sqrRmaxMinusTol;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Check radius and initialize dada members

void G4Orb::Initialize()
{
  const G4double fEpsilon = 2.e-11;  // relative tolerance of fRmax

  // Check radius
  //
  if ( fRmax < 10*kCarTolerance )
  {
    G4Exception("G4Orb::Initialize()", "GeomSolids0002", FatalException,
                "Invalid radius < 10*kCarTolerance.");
  }
  halfRmaxTol = 0.5 * std::max(kCarTolerance, fEpsilon*fRmax);
  G4double rmaxPlusTol  = fRmax + halfRmaxTol;
  G4double rmaxMinusTol = fRmax - halfRmaxTol;
  sqrRmaxPlusTol = rmaxPlusTol*rmaxPlusTol;
  sqrRmaxMinusTol = rmaxMinusTol*rmaxMinusTol;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification

void G4Orb::ComputeDimensions(       G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Orb::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double radius = GetRadius();
  pMin.set(-radius,-radius,-radius);
  pMax.set( radius, radius, radius);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Orb::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Orb::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Check bounding box
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Find bounding envelope and calculate extent
  //
  static const G4int NTHETA = 8;  // number of steps along Theta
  static const G4int NPHI   = 16; // number of steps along Phi
  static const G4double sinHalfTheta = std::sin(halfpi/NTHETA);
  static const G4double cosHalfTheta = std::cos(halfpi/NTHETA);
  static const G4double sinHalfPhi   = std::sin(pi/NPHI);
  static const G4double cosHalfPhi   = std::cos(pi/NPHI);
  static const G4double sinStepTheta = 2.*sinHalfTheta*cosHalfTheta;
  static const G4double cosStepTheta = 1. - 2.*sinHalfTheta*sinHalfTheta;
  static const G4double sinStepPhi   = 2.*sinHalfPhi*cosHalfPhi;
  static const G4double cosStepPhi   = 1. - 2.*sinHalfPhi*sinHalfPhi;

  G4double radius = GetRadius();
  G4double rtheta = radius/cosHalfTheta;
  G4double rphi   = rtheta/cosHalfPhi;

  // set reference circle
  G4TwoVector xy[NPHI];
  G4double sinCurPhi = sinHalfPhi;
  G4double cosCurPhi = cosHalfPhi;
  for (G4int k=0; k<NPHI; ++k)
  {
    xy[k].set(cosCurPhi,sinCurPhi);
    G4double sinTmpPhi = sinCurPhi;
    sinCurPhi = sinCurPhi*cosStepPhi + cosCurPhi*sinStepPhi;
    cosCurPhi = cosCurPhi*cosStepPhi - sinTmpPhi*sinStepPhi;
  }
  
  // set bounding circles
  G4ThreeVectorList circles[NTHETA];
  for (G4int i=0; i<NTHETA; ++i) { circles[i].resize(NPHI); }

  G4double sinCurTheta = sinHalfTheta;
  G4double cosCurTheta = cosHalfTheta;
  for (G4int i=0; i<NTHETA; ++i)
  {
    G4double z = rtheta*cosCurTheta;
    G4double rho = rphi*sinCurTheta;
    for (G4int k=0; k<NPHI; ++k)
    {
      circles[i][k].set(rho*xy[k].x(),rho*xy[k].y(),z);
    }
    G4double sinTmpTheta = sinCurTheta;
    sinCurTheta = sinCurTheta*cosStepTheta + cosCurTheta*sinStepTheta;
    cosCurTheta = cosCurTheta*cosStepTheta - sinTmpTheta*sinStepTheta;
  }

  // set envelope and calculate extent
  std::vector<const G4ThreeVectorList *> polygons;
  polygons.resize(NTHETA);
  for (G4int i=0; i<NTHETA; ++i) { polygons[i] = &circles[i]; }

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Return whether point is inside/outside/on surface

EInside G4Orb::Inside( const G4ThreeVector& p ) const
{
  G4double rr = p.mag2();
  if (rr > sqrRmaxPlusTol) return kOutside;
  return (rr > sqrRmaxMinusTol) ? kSurface : kInside;
}

//////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p

G4ThreeVector G4Orb::SurfaceNormal( const G4ThreeVector& p ) const
{
  return (1/p.mag())*p;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to the surface of the orb from outside
// - return kInfinity if no intersection or
//   intersection distance <= tolerance

G4double G4Orb::DistanceToIn( const G4ThreeVector& p,
                              const G4ThreeVector& v  ) const
{
  // Check if point is on the surface and traveling away
  //
  G4double rr = p.mag2();
  G4double pv = p.dot(v);
  if (rr >= sqrRmaxMinusTol && pv >= 0) return kInfinity;

  // Find intersection
  //
  //    Sphere eqn: x^2 + y^2 + z^2 = R^2
  //
  //    => (px + t*vx)^2 + (py + t*vy)^2 + (pz + t*vz)^2 = R^2
  //    => r^2 + 2t(p.v) + t^2 = R^2
  //    => tmin = -(p.v) - Sqrt((p.v)^2 - (r^2 - R^2))
  //
  G4double D  = pv*pv - rr + fRmax*fRmax;
  if (D < 0) return kInfinity; // no intersection

  G4double sqrtD = std::sqrt(D);
  G4double dist = -pv - sqrtD;

  // Avoid rounding errors due to precision issues seen on 64 bits systems.
  // Split long distances and recompute
  //
  G4double Dmax = 32*fRmax; 
  if (dist > Dmax)
  {
    dist  = dist - 1.e-8*dist - fRmax; // to stay outside after the move
    dist += DistanceToIn(p + dist*v, v);
    return (dist >= kInfinity) ? kInfinity : dist;
  }

  if (sqrtD*2 <= halfRmaxTol) return kInfinity; // touch
  return (dist < halfRmaxTol) ? 0. : dist;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate shortest distance to the boundary from outside
// - Return 0 if point is inside

G4double G4Orb::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double dist = p.mag() - fRmax;
  return (dist > 0) ? dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to the surface of the orb from inside and
// find normal at exit point, if required
// - when leaving the surface, return 0

G4double G4Orb::DistanceToOut( const G4ThreeVector& p,
                               const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool* validNorm,
                                     G4ThreeVector* n ) const
{
  // Check if point is on the surface and traveling away
  //
  G4double rr = p.mag2();
  G4double pv = p.dot(v);
  if (rr >= sqrRmaxMinusTol && pv > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n = p*(1./std::sqrt(rr));
    }
    return 0.;
  }

  // Find intersection
  //
  //    Sphere eqn: x^2 + y^2 + z^2 = R^2
  //
  //    => (px + t*vx)^2 + (py + t*vy)^2 + (pz + t*vz)^2 = R^2
  //    => r^2 + 2t(p.v) + t^2 = R^2
  //    => tmax = -(p.v) + Sqrt((p.v)^2 - (r^2 - R^2))
  //
  G4double D  = pv*pv - rr + fRmax*fRmax;
  G4double tmax = (D <= 0) ? 0. : std::sqrt(D) - pv;
  if (tmax < halfRmaxTol) tmax = 0.;
  if (calcNorm)
  {
    *validNorm = true;
    G4ThreeVector ptmax = p + tmax*v;
    *n = ptmax*(1./ptmax.mag());
  }
  return tmax;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Orb::DistanceToOut( const G4ThreeVector& p ) const
{
#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: " << GetName() << "\n";
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc);
    G4Exception("G4Trap::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
  }
#endif
  G4double dist = fRmax - p.mag();
  return (dist > 0) ? dist : 0.;
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType

G4GeometryType G4Orb::GetEntityType() const
{
  return G4String("G4Orb");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4Orb::Clone() const
{
  return new G4Orb(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Orb::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Orb\n"
     << " Parameters: \n"
     << "    outer radius: " << fRmax/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Orb::GetPointOnSurface() const
{
  return fRmax * G4RandomDirection();
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Orb::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4VisExtent G4Orb::GetExtent() const
{
  return G4VisExtent (-fRmax, fRmax, -fRmax, fRmax, -fRmax, fRmax);
}

G4Polyhedron* G4Orb::CreatePolyhedron () const
{
  return new G4PolyhedronSphere (0., fRmax, 0., 2*pi, 0., pi);
}

#endif
