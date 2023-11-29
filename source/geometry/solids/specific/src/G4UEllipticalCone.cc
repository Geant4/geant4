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
// Implementation for G4UEllipticalCone wrapper class
//
// 13-08-2019 Gabriele Cosmo, CERN
// --------------------------------------------------------------------

#include "G4EllipticalCone.hh"
#include "G4UEllipticalCone.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalConstants.hh"
#include "G4BoundingEnvelope.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UEllipticalCone::G4UEllipticalCone(const G4String& pName,
                                           G4double a,
                                           G4double b,
                                           G4double h,
                                           G4double cut )
  : Base_t(pName, a, b, h, cut)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UEllipticalCone::G4UEllipticalCone( __void__& a )
  : Base_t(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UEllipticalCone::~G4UEllipticalCone() = default;

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UEllipticalCone::G4UEllipticalCone(const G4UEllipticalCone& rhs)
  : Base_t(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UEllipticalCone& G4UEllipticalCone::operator = (const G4UEllipticalCone& rhs)
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

G4double G4UEllipticalCone::GetSemiAxisX() const
{
  return Base_t::GetSemiAxisX();
}

G4double G4UEllipticalCone::GetSemiAxisY() const
{
  return Base_t::GetSemiAxisY();
}

G4double G4UEllipticalCone::GetZMax() const
{
  return Base_t::GetZMax();
}

G4double G4UEllipticalCone::GetZTopCut() const
{
  return Base_t::GetZTopCut();
}

G4double G4UEllipticalCone::GetSemiAxisMax () const
{
  return std::max(GetSemiAxisX(),GetSemiAxisY());
}

G4double G4UEllipticalCone::GetSemiAxisMin () const
{
  return std::min(GetSemiAxisX(),GetSemiAxisY());
}

//////////////////////////////////////////////////////////////////////////
//
// Modifiers

void G4UEllipticalCone::SetSemiAxis(G4double x, G4double y, G4double z)
{
  Base_t::SetParameters(x, y, z, GetZTopCut());
}

void G4UEllipticalCone::SetZCut(G4double newzTopCut)
{
  Base_t::SetParameters(GetSemiAxisX(), GetSemiAxisY(), GetZMax(), newzTopCut);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UEllipticalCone::Clone() const
{
  return new G4UEllipticalCone(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UEllipticalCone::BoundingLimits(G4ThreeVector& pMin,
                                       G4ThreeVector& pMax) const
{
  G4double zcut   = GetZTopCut();
  G4double height = GetZMax(); 
  G4double xmax   = GetSemiAxisX()*(height+zcut);
  G4double ymax   = GetSemiAxisY()*(height+zcut);
  pMin.set(-xmax,-ymax,-zcut);
  pMax.set( xmax, ymax, zcut);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UEllipticalCone::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UEllipticalCone::CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin,bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = pMin < pMax;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  static const G4int NSTEPS = 48; // number of steps for whole circle
  static const G4double ang = twopi/NSTEPS;
  static const G4double sinHalf = std::sin(0.5*ang);
  static const G4double cosHalf = std::cos(0.5*ang);
  static const G4double sinStep = 2.*sinHalf*cosHalf;
  static const G4double cosStep = 1. - 2.*sinHalf*sinHalf;
  G4double zcut   = bmax.z();
  G4double height = GetZMax(); 
  G4double sxmin  = GetSemiAxisX()*(height-zcut)/cosHalf;
  G4double symin  = GetSemiAxisY()*(height-zcut)/cosHalf;
  G4double sxmax  = bmax.x()/cosHalf;
  G4double symax  = bmax.y()/cosHalf;

  G4double sinCur = sinHalf;
  G4double cosCur = cosHalf;
  G4ThreeVectorList baseA(NSTEPS),baseB(NSTEPS);
  for (G4int k=0; k<NSTEPS; ++k)
  {
    baseA[k].set(sxmax*cosCur,symax*sinCur,-zcut);
    baseB[k].set(sxmin*cosCur,symin*sinCur, zcut);
    
    G4double sinTmp = sinCur;
    sinCur = sinCur*cosStep + cosCur*sinStep;
    cosCur = cosCur*cosStep - sinTmp*sinStep;
  }

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;
  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UEllipticalCone::CreatePolyhedron() const
{
  return new G4PolyhedronEllipticalCone(GetSemiAxisX(), GetSemiAxisY(),
                                        GetZMax(), GetZTopCut());
}

#endif  // G4GEOM_USE_USOLIDS
