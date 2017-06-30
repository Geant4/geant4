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
// Implementation for G4UTubs wrapper class
// --------------------------------------------------------------------

#include "G4Tubs.hh"
#include "G4UTubs.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4UTubs::G4UTubs( const G4String& pName,
                        G4double pRMin, G4double pRMax,
                        G4double pDz,
                        G4double pSPhi, G4double pDPhi )
  : Base_t(pName, pRMin, pRMax, pDz, pSPhi, pDPhi)
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTubs::G4UTubs( __void__& a )
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UTubs::~G4UTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UTubs::G4UTubs(const G4UTubs& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UTubs& G4UTubs::operator = (const G4UTubs& rhs) 
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

G4double G4UTubs::GetInnerRadius() const
{
  return rmin();
}
G4double G4UTubs::GetOuterRadius() const
{
  return rmax();
}
G4double G4UTubs::GetZHalfLength() const
{
  return z();
}
G4double G4UTubs::GetStartPhiAngle() const
{
  return sphi();
}
G4double G4UTubs::GetDeltaPhiAngle() const
{
  return dphi();
}
G4double G4UTubs::GetSinStartPhi() const
{
  return std::sin(GetStartPhiAngle());
}
G4double G4UTubs::GetCosStartPhi() const
{
  return std::cos(GetStartPhiAngle());
}
G4double G4UTubs::GetSinEndPhi() const
{
  return std::sin(GetStartPhiAngle()+GetDeltaPhiAngle());
}
G4double G4UTubs::GetCosEndPhi() const
{
  return std::cos(GetStartPhiAngle()+GetDeltaPhiAngle());
}

void G4UTubs::SetInnerRadius(G4double newRMin)
{
  SetRMin(newRMin);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetOuterRadius(G4double newRMax)
{
  SetRMax(newRMax);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetZHalfLength(G4double newDz)
{
  SetDz(newDz);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetStartPhiAngle(G4double newSPhi, G4bool)
{
  SetSPhi(newSPhi);
  fRebuildPolyhedron = true;
}
void G4UTubs::SetDeltaPhiAngle(G4double newDPhi)
{
  SetDPhi(newDPhi);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UTubs::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*(G4Tubs*)this,n,pRep) ;
}

/////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UTubs::Clone() const
{
  return new G4UTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTubs::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double rmin = GetInnerRadius();
  G4double rmax = GetOuterRadius();
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
    G4Exception("G4UTubs::BoundingLimits()", "GeomMgt0001",
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
      G4Exception("G4UTubs::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UTubs::CalculateExtent(const EAxis pAxis,
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
  G4double rmin = GetInnerRadius();
  G4double rmax = GetOuterRadius();
  G4double dz   = GetZHalfLength();
  G4double dphi = GetDeltaPhiAngle();

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
  G4double rext    = rmax/cosHalf;

  // bounding envelope for full cylinder consists of two polygons,
  // in other cases it is a sequence of quadrilaterals
  if (rmin == 0 && dphi == twopi)
  {
    G4double sinCur = sinHalf;
    G4double cosCur = cosHalf;

    G4ThreeVectorList baseA(NSTEPS),baseB(NSTEPS);
    for (G4int k=0; k<NSTEPS; ++k)
    {
      baseA[k].set(rext*cosCur,rext*sinCur,-dz);
      baseB[k].set(rext*cosCur,rext*sinCur, dz);

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
    pols[0][0].set(rmin*cosStart,rmin*sinStart, dz);
    pols[0][1].set(rmin*cosStart,rmin*sinStart,-dz);
    pols[0][2].set(rmax*cosStart,rmax*sinStart,-dz);
    pols[0][3].set(rmax*cosStart,rmax*sinStart, dz);
    for (G4int k=1; k<ksteps+1; ++k)
    {
      pols[k][0].set(rmin*cosCur,rmin*sinCur, dz);
      pols[k][1].set(rmin*cosCur,rmin*sinCur,-dz);
      pols[k][2].set(rext*cosCur,rext*sinCur,-dz);
      pols[k][3].set(rext*cosCur,rext*sinCur, dz);

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    pols[ksteps+1][0].set(rmin*cosEnd,rmin*sinEnd, dz);
    pols[ksteps+1][1].set(rmin*cosEnd,rmin*sinEnd,-dz);
    pols[ksteps+1][2].set(rmax*cosEnd,rmax*sinEnd,-dz);
    pols[ksteps+1][3].set(rmax*cosEnd,rmax*sinEnd, dz);

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
//
G4Polyhedron* G4UTubs::CreatePolyhedron() const
{
  return new G4PolyhedronTubs(GetInnerRadius(),
                              GetOuterRadius(),
                              GetZHalfLength(),
                              GetStartPhiAngle(),
                              GetDeltaPhiAngle());
}

#endif  // G4GEOM_USE_USOLIDS
