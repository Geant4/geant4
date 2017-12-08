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
// Implementation for G4USphere wrapper class
// --------------------------------------------------------------------

#include "G4Sphere.hh"
#include "G4USphere.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4USphere::G4USphere( const G4String& pName,
                            G4double pRmin, G4double pRmax,
                            G4double pSPhi, G4double pDPhi,
                            G4double pSTheta, G4double pDTheta )
  : Base_t(pName, pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta)
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4USphere::G4USphere( __void__& a )
  : Base_t(a)
{
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4USphere::~G4USphere()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4USphere::G4USphere(const G4USphere& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4USphere& G4USphere::operator = (const G4USphere& rhs) 
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

G4double G4USphere::GetInnerRadius() const
{
  return Base_t::GetInnerRadius();
}
G4double G4USphere::GetOuterRadius() const
{
  return Base_t::GetOuterRadius();
}
G4double G4USphere::GetStartPhiAngle() const
{
  return Base_t::GetStartPhiAngle();
}
G4double G4USphere::GetDeltaPhiAngle() const
{
  return Base_t::GetDeltaPhiAngle();
}
G4double G4USphere::GetStartThetaAngle() const
{
  return Base_t::GetStartThetaAngle();
}
G4double G4USphere::GetDeltaThetaAngle() const
{
  return Base_t::GetDeltaThetaAngle();
}
G4double G4USphere::GetSinStartPhi() const
{
  return Base_t::GetSinSPhi();
}
G4double G4USphere::GetCosStartPhi() const
{
  return Base_t::GetCosSPhi();
}
G4double G4USphere::GetSinEndPhi() const
{
  return Base_t::GetSinEPhi();
}
G4double G4USphere::GetCosEndPhi() const
{
  return Base_t::GetCosEPhi();
}
G4double G4USphere::GetSinStartTheta() const
{
  return Base_t::GetSinSTheta();
}
G4double G4USphere::GetCosStartTheta() const
{
  return Base_t::GetCosSTheta();
}
G4double G4USphere::GetSinEndTheta() const
{
  return Base_t::GetSinETheta();
}
G4double G4USphere::GetCosEndTheta() const
{
  return Base_t::GetCosETheta();
}

void G4USphere::SetInnerRadius(G4double newRMin)
{
  Base_t::SetInnerRadius(newRMin);
  fRebuildPolyhedron = true;
}
void G4USphere::SetOuterRadius(G4double newRmax)
{
  Base_t::SetOuterRadius(newRmax);
  fRebuildPolyhedron = true;
}
void G4USphere::SetStartPhiAngle(G4double newSphi, G4bool trig)
{
  Base_t::SetStartPhiAngle(newSphi, trig);
  fRebuildPolyhedron = true;
}
void G4USphere::SetDeltaPhiAngle(G4double newDphi)
{
  Base_t::SetDeltaPhiAngle(newDphi);
  fRebuildPolyhedron = true;
}
void G4USphere::SetStartThetaAngle(G4double newSTheta)
{
  Base_t::SetStartThetaAngle(newSTheta);
  fRebuildPolyhedron = true;
}
void G4USphere::SetDeltaThetaAngle(G4double newDTheta)
{
  Base_t::SetDeltaThetaAngle(newDTheta);
  fRebuildPolyhedron = true;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4USphere::ComputeDimensions(      G4VPVParameterisation* p,
                                  const G4int n,
                                  const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Sphere*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4USphere::Clone() const
{
  return new G4USphere(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4USphere::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double rmin = GetInnerRadius();
  G4double rmax = GetOuterRadius();

  // Find bounding box
  //
  if (GetDeltaThetaAngle() >= pi && GetDeltaPhiAngle() >= twopi)
  {
    pMin.set(-rmax,-rmax,-rmax);
    pMax.set( rmax, rmax, rmax);
  }
  else
  {
    G4double sinStart = GetSinStartTheta();
    G4double cosStart = GetCosStartTheta();
    G4double sinEnd   = GetSinEndTheta();
    G4double cosEnd   = GetCosEndTheta();

    G4double stheta = GetStartThetaAngle();
    G4double etheta = stheta + GetDeltaThetaAngle();
    G4double rhomin = rmin*std::min(sinStart,sinEnd);
    G4double rhomax = rmax;
    if (stheta > halfpi) rhomax = rmax*sinStart;
    if (etheta < halfpi) rhomax = rmax*sinEnd;

    G4TwoVector xymin,xymax;
    G4GeomTools::DiskExtent(rhomin,rhomax,
                            GetSinStartPhi(),GetCosStartPhi(),
                            GetSinEndPhi(),GetCosEndPhi(),
                            xymin,xymax);

    G4double zmin = std::min(rmin*cosEnd,rmax*cosEnd);
    G4double zmax = std::max(rmin*cosStart,rmax*cosStart);
    pMin.set(xymin.x(),xymin.y(),zmin);
    pMax.set(xymax.x(),xymax.y(),zmax);
  }

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4USphere::BoundingLimits()", "GeomMgt0001",
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
      G4Exception("G4USphere::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4USphere::CalculateExtent(const EAxis pAxis,
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
// Create polyhedron for visualization

G4Polyhedron* G4USphere::CreatePolyhedron() const
{
  return new G4PolyhedronSphere(GetInnerRadius(),
                                GetOuterRadius(),
                                GetStartPhiAngle(),
                                GetDeltaPhiAngle(),
                                GetStartThetaAngle(),
                                GetDeltaThetaAngle());
}

#endif  // G4GEOM_USE_USOLIDS
