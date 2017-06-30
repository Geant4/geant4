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
// Implementation for G4UParaboloid wrapper class
//
// 19-08-2015 Guilherme Lima, FNAL
//
// --------------------------------------------------------------------

#include "G4Paraboloid.hh"
#include "G4UParaboloid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4BoundingEnvelope.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UParaboloid::G4UParaboloid(const G4String& pName,
                                   G4double dz,
                                   G4double rlo,
                                   G4double rhi )
  : Base_t(pName, rlo, rhi, dz)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UParaboloid::G4UParaboloid( __void__& a )
  : Base_t(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UParaboloid::~G4UParaboloid() { }

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UParaboloid::G4UParaboloid(const G4UParaboloid& rhs)
  : Base_t(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UParaboloid& G4UParaboloid::operator = (const G4UParaboloid& rhs)
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
// Accessors

G4double G4UParaboloid::GetZHalfLength() const
{
  return GetDz();
}

G4double G4UParaboloid::GetRadiusMinusZ() const
{
  return GetRlo();
}

G4double G4UParaboloid::GetRadiusPlusZ() const
{
  return GetRhi();
}

//////////////////////////////////////////////////////////////////////////
//
// Modifiers

void G4UParaboloid::SetZHalfLength(G4double dz)
{
  SetDz(dz);
}

void G4UParaboloid::SetRadiusMinusZ(G4double r1)
{
  SetRlo(r1);
}

void G4UParaboloid::SetRadiusPlusZ(G4double r2)
{
  SetRhi(r2);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UParaboloid::Clone() const
{
  return new G4UParaboloid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UParaboloid::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double r2 = GetRadiusPlusZ();
  G4double dz = GetZHalfLength();
  pMin.set(-r2,-r2,-dz);
  pMax.set( r2, r2, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UParaboloid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  // Check consistency of bounding boxes
  //
  if (checkBBox)
  {
    U3Vector vmin, vmax;
    Extent(vmin,vmax);
    if (std::abs(pMin.x()-vmin.x()) > kCarTolerance ||
        std::abs(pMin.y()-vmin.y()) > kCarTolerance ||
        std::abs(pMin.z()-vmin.z()) > kCarTolerance ||
        std::abs(pMax.x()-vmax.x()) > kCarTolerance ||
        std::abs(pMax.y()-vmax.y()) > kCarTolerance ||
        std::abs(pMax.z()-vmax.z()) > kCarTolerance)
    {
      std::ostringstream message;
      message << "Inconsistency in bounding boxes for solid: "
              << GetName() << " !"
              << "\nBBox min: wrapper = " << pMin << " solid = " << vmin
              << "\nBBox max: wrapper = " << pMax << " solid = " << vmax;
      G4Exception("G4UParaboloid::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UParaboloid::CalculateExtent(const EAxis pAxis,
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

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UParaboloid::CreatePolyhedron() const
{
  return new G4PolyhedronParaboloid(GetRadiusMinusZ(),
                                    GetRadiusPlusZ(),
                                    GetZHalfLength(), 0., twopi);
}

#endif  // G4GEOM_USE_USOLIDS
