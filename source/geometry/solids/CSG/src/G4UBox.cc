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
// 
// Implementation for G4UBox wrapper class
// --------------------------------------------------------------------

#include "G4Box.hh"
#include "G4UBox.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UBox::G4UBox(const G4String& pName,
                     G4double pX,
                     G4double pY,
                     G4double pZ)
  : Base_t(pName, pX, pY, pZ)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UBox::G4UBox( __void__& a )
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UBox::~G4UBox()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UBox::G4UBox(const G4UBox& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UBox& G4UBox::operator = (const G4UBox& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   Base_t::operator=(rhs);

   return *this;
}

////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers

G4double G4UBox::GetXHalfLength() const
{
  return x();
}
G4double G4UBox::GetYHalfLength() const
{
  return y();
}
G4double G4UBox::GetZHalfLength() const
{
  return z();
}

void G4UBox::SetXHalfLength(G4double dx)
{
  SetX(dx);
  fRebuildPolyhedron = true;
}
void G4UBox::SetYHalfLength(G4double dy)
{
  SetY(dy);
  fRebuildPolyhedron = true;
}
void G4UBox::SetZHalfLength(G4double dz)
{
  SetZ(dz);
  fRebuildPolyhedron = true;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UBox::ComputeDimensions(G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Box*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UBox::Clone() const
{
  return new G4UBox(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// StreamInfo

std::ostream& G4UBox::StreamInfo(std::ostream &os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "     *** Dump for solid - " << GetEntityType() << " ***\n"
     << "     ===================================================\n"
     << " Solid type: Box\n"
     << " Parameters: \n"
     << "     half-dimensions in mm: x,y,z: " << dimensions() <<"\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UBox::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double dx = GetXHalfLength();
  G4double dy = GetYHalfLength();
  G4double dz = GetZHalfLength();
  pMin.set(-dx,-dy,-dz);
  pMax.set( dx, dy, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UBox::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UBox::CalculateExtent(const EAxis pAxis,
                        const G4VoxelLimits& pVoxelLimit,
                        const G4AffineTransform& pTransform,
                              G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box limits
  BoundingLimits(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}


//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization

G4Polyhedron* G4UBox::CreatePolyhedron() const
{
  return new G4PolyhedronBox(GetXHalfLength(),
                             GetYHalfLength(),
                             GetZHalfLength());
}

#endif  // G4GEOM_USE_USOLIDS
