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
// Implementation for G4UTrd wrapper class
// --------------------------------------------------------------------

#include "G4Trd.hh"
#include "G4UTrd.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths
//
G4UTrd::G4UTrd(const G4String& pName,
                     G4double pdx1,  G4double pdx2,
                     G4double pdy1,  G4double pdy2,
                     G4double pdz)
  : Base_t(pName, pdx1, pdx2, pdy1, pdy2, pdz)
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTrd::G4UTrd( __void__& a )
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTrd::~G4UTrd()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTrd::G4UTrd(const G4UTrd& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTrd& G4UTrd::operator = (const G4UTrd& rhs) 
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

G4double G4UTrd::GetXHalfLength1() const
{
  return dx1();
}
G4double G4UTrd::GetXHalfLength2() const
{
  return dx2();
}
G4double G4UTrd::GetYHalfLength1() const
{
  return dy1();
}
G4double G4UTrd::GetYHalfLength2() const
{
  return dy2();
}
G4double G4UTrd::GetZHalfLength()  const
{
  return dz();
}

void G4UTrd::SetXHalfLength1(G4double val)
{
  Base_t::SetXHalfLength1(val);
  fRebuildPolyhedron = true;
}
void G4UTrd::SetXHalfLength2(G4double val)
{
  Base_t::SetXHalfLength2(val);
  fRebuildPolyhedron = true;
}
void G4UTrd::SetYHalfLength1(G4double val)
{
  Base_t::SetYHalfLength1(val);
  fRebuildPolyhedron = true;
}
void G4UTrd::SetYHalfLength2(G4double val)
{
  Base_t::SetYHalfLength2(val);
  fRebuildPolyhedron = true;
}
void G4UTrd::SetZHalfLength(G4double val)
{
  Base_t::SetZHalfLength(val);
  fRebuildPolyhedron = true;
}
void G4UTrd::SetAllParameters(G4double pdx1, G4double pdx2,
                              G4double pdy1, G4double pdy2, G4double pdz)
{
  Base_t::SetAllParameters(pdx1, pdx2, pdy1, pdy2, pdz);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4UTrd::ComputeDimensions(      G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Trd*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UTrd::Clone() const
{
  return new G4UTrd(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTrd::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double dx1 = GetXHalfLength1();
  G4double dx2 = GetXHalfLength2();
  G4double dy1 = GetYHalfLength1();
  G4double dy2 = GetYHalfLength2();
  G4double dz  = GetZHalfLength();

  G4double xmax = std::max(dx1,dx2);
  G4double ymax = std::max(dy1,dy2);
  pMin.set(-xmax,-ymax,-dz);
  pMax.set( xmax, ymax, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UTrd::BoundingLimits()", "GeomMgt0001",
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
      G4Exception("G4UTrd::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UTrd::CalculateExtent(const EAxis pAxis,
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
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  G4double dx1 = GetXHalfLength1();
  G4double dx2 = GetXHalfLength2();
  G4double dy1 = GetYHalfLength1();
  G4double dy2 = GetYHalfLength2();
  G4double dz  = GetZHalfLength();

  G4ThreeVectorList baseA(4), baseB(4);
  baseA[0].set(-dx1,-dy1,-dz);
  baseA[1].set( dx1,-dy1,-dz);
  baseA[2].set( dx1, dy1,-dz);
  baseA[3].set(-dx1, dy1,-dz);
  baseB[0].set(-dx2,-dy2, dz);
  baseB[1].set( dx2,-dy2, dz);
  baseB[2].set( dx2, dy2, dz);
  baseB[3].set(-dx2, dy2, dz);

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization
//
G4Polyhedron* G4UTrd::CreatePolyhedron() const
{
  return new G4PolyhedronTrd2(GetXHalfLength1(),
                              GetXHalfLength2(),
                              GetYHalfLength1(),
                              GetYHalfLength2(),
                              GetZHalfLength());
}

#endif  // G4GEOM_USE_USOLIDS
