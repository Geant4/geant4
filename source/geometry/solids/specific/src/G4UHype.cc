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
// Implementation for G4UHype wrapper class
//
// 16-10-2017 G.Cosmo, CERN
//
// --------------------------------------------------------------------

#include "G4Hype.hh"
#if 0
#include "G4UHype.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor

G4UHype::G4UHype(const G4String& pName,
                       G4double  newInnerRadius,
                       G4double  newOuterRadius,
                       G4double  newInnerStereo,
                       G4double  newOuterStereo,
                       G4double  newHalfLenZ)
  : Base_t(pName, newInnerRadius, newOuterRadius,
                  newInnerStereo, newOuterStereo, newHalfLenZ)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UHype::G4UHype( __void__& a )
  : Base_t(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UHype::~G4UHype() { }

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UHype::G4UHype(const G4UHype& rhs)
  : Base_t(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UHype& G4UHype::operator = (const G4UHype& rhs)
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

G4double G4UHype::GetInnerRadius () const
{
  return GetRmin();
}

G4double G4UHype::GetOuterRadius () const
{
  return GetRmax();
}

G4double G4UHype::GetZHalfLength () const
{
  return GetDz();
}

G4double G4UHype::GetInnerStereo () const
{
  return GetStIn();
}

G4double G4UHype::GetOuterStereo () const
{
  return GetStOut();
}

//////////////////////////////////////////////////////////////////////////
//
// Modifiers

void G4UHype::SetInnerRadius (G4double newIRad)
{
  SetParameters(newIRad, GetRmax(), GetStIn(), GetStOut(), GetDz());
  fRebuildPolyhedron = true;
}

void G4UHype::SetOuterRadius (G4double newORad)
{
  SetParameters(GetRmin(), newORad, GetStIn(), GetStOut(), GetDz());
  fRebuildPolyhedron = true;
}

void G4UHype::SetZHalfLength (G4double newHLZ)
{
  SetParameters(GetRmin(), GetRmax(), GetStIn(), GetStOut(), newHLZ);
  fRebuildPolyhedron = true;
}

void G4UHype::SetInnerStereo (G4double newISte)
{
  SetParameters(GetRmin(), GetRmax(), newISte, GetStOut(), GetDz());
  fRebuildPolyhedron = true;
}

void G4UHype::SetOuterStereo (G4double newOSte)
{
  SetParameters(GetRmin(), GetRmax(), GetStIn(), newOSte, GetDz());
  fRebuildPolyhedron = true;
}


////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UHype::ComputeDimensions(G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Hype*)this,n,pRep);
}


//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UHype::Clone() const
{
  return new G4UHype(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UHype::BoundingLimits(G4ThreeVector& pMin,
                             G4ThreeVector& pMax) const
{
  G4double endORadius = GetEndInnerRadius();
  pMin.set(-endORadius,-endORadius,-GetDz());
  pMax.set( endORadius, endORadius, GetDz());

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UHype::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UHype::CalculateExtent(const EAxis pAxis,
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
G4Polyhedron* G4UHype::CreatePolyhedron() const
{
  return new G4PolyhedronHype(GetRmin(), GetRmax(),
                              GetTIn2(), GetTOut2(), GetDz());
}

#endif  // G4GEOM_USE_USOLIDS

#endif
