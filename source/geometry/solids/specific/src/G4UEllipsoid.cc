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
// Implementation for G4UEllipsoid wrapper class
//
// 13-08-2019 Gabriele Cosmo, CERN
// --------------------------------------------------------------------

#include "G4Ellipsoid.hh"
#include "G4UEllipsoid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4BoundingEnvelope.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UEllipsoid::G4UEllipsoid(const G4String& pName,
                                           G4double dx,
                                           G4double dy,
                                           G4double dz,
                                           G4double bcut,
                                           G4double tcut )
  : Base_t(pName, dx, dy, dz, bcut, tcut)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UEllipsoid::G4UEllipsoid( __void__& a )
  : Base_t(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UEllipsoid::~G4UEllipsoid() = default;

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UEllipsoid::G4UEllipsoid(const G4UEllipsoid& rhs)
  : Base_t(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UEllipsoid& G4UEllipsoid::operator = (const G4UEllipsoid& rhs)
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

G4double G4UEllipsoid::GetDx() const
{
  return Base_t::GetDx();
}

G4double G4UEllipsoid::GetDy() const
{
  return Base_t::GetDy();
}

G4double G4UEllipsoid::GetDz() const
{
  return Base_t::GetDz();
}

G4double G4UEllipsoid::GetSemiAxisMax (G4int i) const
{
  return (i==0) ? GetDx()
       : (i==1) ? GetDy()
       : GetDz();
}

G4double G4UEllipsoid::GetZBottomCut() const
{
  return Base_t::GetZBottomCut();
}

G4double G4UEllipsoid::GetZTopCut() const
{
  return Base_t::GetZTopCut();
}

//////////////////////////////////////////////////////////////////////////
//
// Modifiers

void G4UEllipsoid::SetSemiAxis (G4double x, G4double y, G4double z)
{
  Base_t::SetSemiAxes(x, y, z);
}

void G4UEllipsoid::SetZCuts (G4double newzBottomCut, G4double newzTopCut)
{
  Base_t::SetZCuts(newzBottomCut, newzTopCut);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UEllipsoid::Clone() const
{
  return new G4UEllipsoid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UEllipsoid::BoundingLimits(G4ThreeVector& pMin,
                                  G4ThreeVector& pMax) const
{
  G4double dx = GetDx();
  G4double dy = GetDy();
  G4double dz = GetDz();
  G4double zmin = std::max(-dz,GetZBottomCut());
  G4double zmax = std::min( dz,GetZTopCut());
  pMin.set(-dx,-dy,zmin);
  pMax.set( dx, dy,zmax);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UEllipsoid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UEllipsoid::CalculateExtent(const EAxis pAxis,
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

G4Polyhedron* G4UEllipsoid::CreatePolyhedron() const
{
  return new G4PolyhedronEllipsoid(GetDx(), GetDy(), GetDz(),
                                   GetZBottomCut(), GetZTopCut());
}

#endif  // G4GEOM_USE_USOLIDS
