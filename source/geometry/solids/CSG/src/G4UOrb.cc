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
// $Id:$
//
// 
// Implementation for G4UOrb wrapper class
// --------------------------------------------------------------------

#include "G4Orb.hh"
#include "G4UOrb.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4TwoVector.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"
#include "G4PhysicalConstants.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// constructor - check positive radius
//             

G4UOrb::G4UOrb( const G4String& pName, G4double pRmax )
  : Base_t(pName, pRmax)
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UOrb::G4UOrb( __void__& a )
  : Base_t(a)
{
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4UOrb::~G4UOrb()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UOrb::G4UOrb(const G4UOrb& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UOrb& G4UOrb::operator = (const G4UOrb& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   Base_t::operator=(rhs);

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers

G4double G4UOrb::GetRadius() const
{
  return Base_t::GetRadius();
}

void G4UOrb::SetRadius(G4double newRmax)
{
  Base_t::SetRadius(newRmax);
  fRebuildPolyhedron = true;
}

G4double G4UOrb::GetRadialTolerance() const
{
  return Base_t::GetRadialTolerance();
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UOrb::ComputeDimensions(       G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*(G4Orb*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UOrb::Clone() const
{
  return new G4UOrb(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UOrb::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
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
    G4Exception("G4UOrb::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UOrb::CalculateExtent(const EAxis pAxis,
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
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
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
  for (G4int i=0; i<NTHETA; ++i) circles[i].resize(NPHI);

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
  for (G4int i=0; i<NTHETA; ++i) polygons[i] = &circles[i];

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization

G4Polyhedron* G4UOrb::CreatePolyhedron() const
{
  return new G4PolyhedronSphere(0., GetRadius(), 0., twopi, 0., pi);
}

#endif  // G4GEOM_USE_USOLIDS
