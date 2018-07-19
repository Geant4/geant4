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
// Implementation for G4UCons wrapper class
// --------------------------------------------------------------------

#include "G4Cons.hh"
#include "G4UCons.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

G4UCons::G4UCons( const G4String& pName,
                        G4double  pRmin1, G4double pRmax1,
                        G4double  pRmin2, G4double pRmax2,
                        G4double pDz,
                        G4double pSPhi, G4double pDPhi)
  : Base_t(pName, pRmin1, pRmax1, pRmin2, pRmax2, pDz, pSPhi, pDPhi)
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UCons::G4UCons( __void__& a )
  : Base_t(a)
{
}

///////////////////////////////////////////////////////////////////////
//
// Destructor

G4UCons::~G4UCons()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UCons::G4UCons(const G4UCons& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UCons& G4UCons::operator = (const G4UCons& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   Base_t::operator=(rhs);

   return *this;
}

/////////////////////////////////////////////////////////////////////////
//
// Accessors and modifiers

G4double G4UCons::GetInnerRadiusMinusZ() const
{
  return GetRmin1();
}
G4double G4UCons::GetOuterRadiusMinusZ() const
{
  return GetRmax1();
}
G4double G4UCons::GetInnerRadiusPlusZ() const
{
  return GetRmin2();
}
G4double G4UCons::GetOuterRadiusPlusZ() const
{
  return GetRmax2();
}
G4double G4UCons::GetZHalfLength() const
{
  return GetDz();
}
G4double G4UCons::GetStartPhiAngle() const
{
  return GetSPhi();
}
G4double G4UCons::GetDeltaPhiAngle() const
{
  return GetDPhi();
}
G4double G4UCons::GetSinStartPhi() const
{
  G4double phi = GetStartPhiAngle();
  return std::sin(phi);
}
G4double G4UCons::GetCosStartPhi() const
{
  G4double phi = GetStartPhiAngle();
  return std::cos(phi);
}
G4double G4UCons::GetSinEndPhi() const
{
  G4double phi = GetStartPhiAngle() + GetDeltaPhiAngle();
  return std::sin(phi);
}
G4double G4UCons::GetCosEndPhi() const
{
  G4double phi = GetStartPhiAngle() + GetDeltaPhiAngle();
  return std::cos(phi);
}
  
void G4UCons::SetInnerRadiusMinusZ(G4double Rmin1)
{
  SetRmin1(Rmin1);
  fRebuildPolyhedron = true;
}
void G4UCons::SetOuterRadiusMinusZ(G4double Rmax1)
{
  SetRmax1(Rmax1);
  fRebuildPolyhedron = true;
}
void G4UCons::SetInnerRadiusPlusZ(G4double Rmin2)
{
  SetRmin2(Rmin2);
  fRebuildPolyhedron = true;
}
void G4UCons::SetOuterRadiusPlusZ(G4double Rmax2)
{
  SetRmax2(Rmax2);
  fRebuildPolyhedron = true;
}
void G4UCons::SetZHalfLength(G4double newDz)
{
  SetDz(newDz);
  fRebuildPolyhedron = true;
}
void G4UCons::SetStartPhiAngle(G4double newSPhi, G4bool)
{
  SetSPhi(newSPhi);
  fRebuildPolyhedron = true;
}
void G4UCons::SetDeltaPhiAngle(G4double newDPhi)
{
  SetDPhi(newDPhi);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UCons::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int                  n,
                                const G4VPhysicalVolume*     pRep    )
{
  p->ComputeDimensions(*(G4Cons*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UCons::Clone() const
{
  return new G4UCons(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UCons::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double rmin = std::min(GetInnerRadiusMinusZ(),GetInnerRadiusPlusZ());
  G4double rmax = std::max(GetOuterRadiusMinusZ(),GetOuterRadiusPlusZ());
  G4double dz   = GetZHalfLength();

  // Find bounding box
  //
  if (GetDeltaPhiAngle() < twopi)
  {
    G4TwoVector vmin,vmax;
    G4GeomTools::DiskExtent(rmin,rmax,
                            GetSinStartPhi(),GetCosStartPhi(),
                            GetSinEndPhi(),GetCosEndPhi(),
                            vmin,vmax);
    pMin.set(vmin.x(),vmin.y(),-dz);
    pMax.set(vmax.x(),vmax.y(), dz);
  }
  else
  {
    pMin.set(-rmax,-rmax,-dz);
    pMax.set( rmax, rmax, dz);
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
    G4Exception("G4UCons::BoundingLimits()", "GeomMgt0001",
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
      G4Exception("G4UCons::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UCons::CalculateExtent(const EAxis pAxis,
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

  // Get parameters of the solid
  G4double rmin1 = GetInnerRadiusMinusZ();
  G4double rmax1 = GetOuterRadiusMinusZ();
  G4double rmin2 = GetInnerRadiusPlusZ();
  G4double rmax2 = GetOuterRadiusPlusZ();
  G4double dz    = GetZHalfLength();
  G4double dphi  = GetDeltaPhiAngle();

  // Find bounding envelope and calculate extent
  //
  const G4int NSTEPS = 24;            // number of steps for whole circle
  G4double astep  = twopi/NSTEPS;     // max angle for one step
  G4int    ksteps = (dphi <= astep) ? 1 : (G4int)((dphi-deg)/astep) + 1;
  G4double ang    = dphi/ksteps;

  G4double sinHalf = std::sin(0.5*ang);
  G4double cosHalf = std::cos(0.5*ang);
  G4double sinStep = 2.*sinHalf*cosHalf;
  G4double cosStep = 1. - 2.*sinHalf*sinHalf;
  G4double rext1   = rmax1/cosHalf;
  G4double rext2   = rmax2/cosHalf;

  // bounding envelope for full cone without hole consists of two polygons,
  // in other cases it is a sequence of quadrilaterals
  if (rmin1 == 0 && rmin2 == 0 && dphi == twopi)
  {
    G4double sinCur = sinHalf;
    G4double cosCur = cosHalf;

    G4ThreeVectorList baseA(NSTEPS),baseB(NSTEPS);
    for (G4int k=0; k<NSTEPS; ++k)
    {
      baseA[k].set(rext1*cosCur,rext1*sinCur,-dz);
      baseB[k].set(rext2*cosCur,rext2*sinCur, dz);

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    std::vector<const G4ThreeVectorList *> polygons(2);
    polygons[0] = &baseA;
    polygons[1] = &baseB;
    G4BoundingEnvelope benv(bmin,bmax,polygons);
    exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }
  else
  {
    G4double sinStart = GetSinStartPhi();
    G4double cosStart = GetCosStartPhi();
    G4double sinEnd   = GetSinEndPhi();
    G4double cosEnd   = GetCosEndPhi();
    G4double sinCur   = sinStart*cosHalf + cosStart*sinHalf;
    G4double cosCur   = cosStart*cosHalf - sinStart*sinHalf;

    // set quadrilaterals
    G4ThreeVectorList pols[NSTEPS+2];
    for (G4int k=0; k<ksteps+2; ++k) pols[k].resize(4);
    pols[0][0].set(rmin2*cosStart,rmin2*sinStart, dz);
    pols[0][1].set(rmin1*cosStart,rmin1*sinStart,-dz);
    pols[0][2].set(rmax1*cosStart,rmax1*sinStart,-dz);
    pols[0][3].set(rmax2*cosStart,rmax2*sinStart, dz);
    for (G4int k=1; k<ksteps+1; ++k)
    {
      pols[k][0].set(rmin2*cosCur,rmin2*sinCur, dz);
      pols[k][1].set(rmin1*cosCur,rmin1*sinCur,-dz);
      pols[k][2].set(rext1*cosCur,rext1*sinCur,-dz);
      pols[k][3].set(rext2*cosCur,rext2*sinCur, dz);

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    pols[ksteps+1][0].set(rmin2*cosEnd,rmin2*sinEnd, dz);
    pols[ksteps+1][1].set(rmin1*cosEnd,rmin1*sinEnd,-dz);
    pols[ksteps+1][2].set(rmax1*cosEnd,rmax1*sinEnd,-dz);
    pols[ksteps+1][3].set(rmax2*cosEnd,rmax2*sinEnd, dz);

    // set envelope and calculate extent
    std::vector<const G4ThreeVectorList *> polygons;
    polygons.resize(ksteps+2);
    for (G4int k=0; k<ksteps+2; ++k) polygons[k] = &pols[k];
    G4BoundingEnvelope benv(bmin,bmax,polygons);
    exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization

G4Polyhedron* G4UCons::CreatePolyhedron() const
{
  return new G4PolyhedronCons(GetInnerRadiusMinusZ(),
                              GetOuterRadiusMinusZ(),
                              GetInnerRadiusPlusZ(),
                              GetOuterRadiusPlusZ(),
                              GetZHalfLength(),
                              GetStartPhiAngle(),
                              GetDeltaPhiAngle());
}

#endif  // G4GEOM_USE_USOLIDS
