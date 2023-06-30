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
// Implementation for G4UEllipticalTube wrapper class
//
// 13-08-2019 Gabriele Cosmo, CERN
// --------------------------------------------------------------------

#include "G4EllipticalTube.hh"
#include "G4UEllipticalTube.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4BoundingEnvelope.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UEllipticalTube::G4UEllipticalTube(const G4String& pName,
                                           G4double dx,
                                           G4double dy,
                                           G4double dz )
  : Base_t(pName, dx, dy, dz)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UEllipticalTube::G4UEllipticalTube( __void__& a )
  : Base_t(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UEllipticalTube::~G4UEllipticalTube() = default;

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UEllipticalTube::G4UEllipticalTube(const G4UEllipticalTube& rhs)
  : Base_t(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UEllipticalTube& G4UEllipticalTube::operator = (const G4UEllipticalTube& rhs)
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

G4double G4UEllipticalTube::GetDx() const
{
  return Base_t::GetDx();
}

G4double G4UEllipticalTube::GetDy() const
{
  return Base_t::GetDy();
}

G4double G4UEllipticalTube::GetDz() const
{
  return Base_t::GetDz();
}

//////////////////////////////////////////////////////////////////////////
//
// Modifiers

void G4UEllipticalTube::SetDx(G4double dx)
{
  Base_t::SetDx(dx);
}

void G4UEllipticalTube::SetDy(G4double dy)
{
  Base_t::SetDy(dy);
}

void G4UEllipticalTube::SetDz(G4double dz)
{
  Base_t::SetDz(dz);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UEllipticalTube::Clone() const
{
  return new G4UEllipticalTube(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UEllipticalTube::BoundingLimits(G4ThreeVector& pMin,
                                       G4ThreeVector& pMax) const
{
  G4double dx = GetDx();
  G4double dy = GetDy();
  G4double dz = GetDz();

  pMin.set(-dx,-dy,-dz);
  pMax.set( dx, dy, dz);
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UEllipticalTube::CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double& pMin, G4double& pMax) const
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
    return exist = pMin < pMax;
  }

  G4double dx = GetDx();
  G4double dy = GetDy();
  G4double dz = GetDz();

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
  G4ThreeVectorList baseA(NSTEPS), baseB(NSTEPS);
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

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UEllipticalTube::CreatePolyhedron() const
{
  // create cylinder with radius=1...
  //
  G4Polyhedron* eTube = new G4PolyhedronTube(0., 1., GetDz());

  // apply non-uniform scaling...
  //
  eTube->Transform(G4Scale3D(GetDx(), GetDy(), 1.));
  return eTube;
}

#endif  // G4GEOM_USE_USOLIDS
