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
// Implementation for G4UTorus wrapper class
//
// 19-08-2015 Guilherme Lima, FNAL
//
// --------------------------------------------------------------------

#include "G4Torus.hh"
#include "G4UTorus.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4TwoVector.hh"
#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths


G4UTorus::G4UTorus(const G4String& pName,
                         G4double rmin, G4double rmax, G4double rtor,
                         G4double sphi, G4double dphi)
  : G4USolid(pName, new UTorus(pName, rmin, rmax, rtor, sphi, dphi))
{ }

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UTorus::G4UTorus( __void__& a )
  : G4USolid(a)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UTorus::~G4UTorus() { }

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UTorus::G4UTorus(const G4UTorus& rhs)
  : G4USolid(rhs)
{ }

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UTorus& G4UTorus::operator = (const G4UTorus& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4USolid::operator=(rhs);

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers

G4double G4UTorus::GetRmin() const
{
  return GetShape()->GetRmin();
}

G4double G4UTorus::GetRmax() const
{
  return GetShape()->GetRmax();
}

G4double G4UTorus::GetRtor() const
{
  return GetShape()->GetRtor();
}

G4double G4UTorus::GetSPhi() const
{
  return GetShape()->GetSPhi();
}

G4double G4UTorus::GetDPhi() const
{
  return GetShape()->GetDPhi();
}

G4double G4UTorus::GetSinStartPhi() const
{
  G4double phi = GetShape()->GetSPhi();
  return std::sin(phi);
}

G4double G4UTorus::GetCosStartPhi() const
{
  G4double phi = GetShape()->GetSPhi();
  return std::cos(phi);
}

G4double G4UTorus::GetSinEndPhi() const
{
  G4double phi = GetShape()->GetSPhi() +
                 GetShape()->GetDPhi();
  return std::sin(phi);
}

G4double G4UTorus::GetCosEndPhi() const
{
  G4double phi = GetShape()->GetSPhi() +
                 GetShape()->GetDPhi();
  return std::cos(phi);
}

void G4UTorus::SetRmin(G4double arg)
{
  GetShape()->SetRmin(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetRmax(G4double arg)
{
  GetShape()->SetRmax(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetRtor(G4double arg)
{
  GetShape()->SetRtor(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetSPhi(G4double arg)
{
  GetShape()->SetSPhi(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetDPhi(G4double arg)
{
  GetShape()->SetDPhi(arg);
  fRebuildPolyhedron = true;
}

void G4UTorus::SetAllParameters(G4double arg1, G4double arg2,
                        G4double arg3, G4double arg4, G4double arg5)
{
  GetShape()->SetRmin(arg1);
  GetShape()->SetRmax(arg2);
  GetShape()->SetRtor(arg3);
  GetShape()->SetSPhi(arg4);
  GetShape()->SetDPhi(arg5);
  fRebuildPolyhedron = true;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4UTorus::ComputeDimensions(G4VPVParameterisation* p,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Torus*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UTorus::Clone() const
{
  return new G4UTorus(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTorus::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double rmax = GetRmax();
  G4double rtor = GetRtor();
  G4double rint = rtor - rmax;
  G4double rext = rtor + rmax;
  G4double dz   = rmax;

  // Find bounding box
  //
  if (GetDPhi() >= twopi)
  {
    pMin.set(-rext,-rext,-dz);
    pMax.set( rext, rext, dz);
  }
  else
  {
    G4TwoVector vmin,vmax;
    G4GeomTools::DiskExtent(rint,rext,
                            GetSinStartPhi(),GetCosStartPhi(),
                            GetSinEndPhi(),GetCosEndPhi(),
                            vmin,vmax);
    pMin.set(vmin.x(),vmin.y(),-dz);
    pMax.set(vmax.x(),vmax.y(), dz);
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
    G4Exception("G4UTorus::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  // Check consistency of bounding boxes
  //
  if (checkBBox)
  {
    UVector3 vmin, vmax;
    GetShape()->Extent(vmin,vmax);
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
      G4Exception("G4UTorus::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UTorus::CalculateExtent(const EAxis pAxis,
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
  G4double rmin = GetRmin();
  G4double rmax = GetRmax();
  G4double rtor = GetRtor();
  G4double dphi = GetDPhi();
  G4double sinStart = GetSinStartPhi();
  G4double cosStart = GetCosStartPhi();
  G4double sinEnd   = GetSinEndPhi();
  G4double cosEnd   = GetCosEndPhi();
  G4double rint = rtor - rmax;
  G4double rext = rtor + rmax;

  // Find bounding envelope and calculate extent
  //
  static const G4int NPHI  = 24; // number of steps for whole torus
  static const G4int NDISK = 16; // number of steps for disk
  static const G4double sinHalfDisk = std::sin(pi/NDISK);
  static const G4double cosHalfDisk = std::cos(pi/NDISK);
  static const G4double sinStepDisk = 2.*sinHalfDisk*cosHalfDisk;
  static const G4double cosStepDisk = 1. - 2.*sinHalfDisk*sinHalfDisk;

  G4double astep = (360/NPHI)*deg; // max angle for one slice in phi
  G4int    kphi  = (dphi <= astep) ? 1 : (G4int)((dphi-deg)/astep) + 1;
  G4double ang   = dphi/kphi;

  G4double sinHalf = std::sin(0.5*ang);
  G4double cosHalf = std::cos(0.5*ang);
  G4double sinStep = 2.*sinHalf*cosHalf;
  G4double cosStep = 1. - 2.*sinHalf*sinHalf;

  // define vectors for bounding envelope
  G4ThreeVectorList pols[NDISK+1];
  for (G4int k=0; k<NDISK+1; ++k) pols[k].resize(4);

  std::vector<const G4ThreeVectorList *> polygons;
  polygons.resize(NDISK+1);
  for (G4int k=0; k<NDISK+1; ++k) polygons[k] = &pols[k];

  // set internal and external reference circles
  G4TwoVector rzmin[NDISK];
  G4TwoVector rzmax[NDISK];

  if ((rtor-rmin*sinHalfDisk)/cosHalf > (rtor+rmin*sinHalfDisk)) rmin = 0;
  rmax /= cosHalfDisk;
  G4double sinCurDisk = sinHalfDisk;
  G4double cosCurDisk = cosHalfDisk;
  for (G4int k=0; k<NDISK; ++k)
  {
    G4double rmincur = rtor + rmin*cosCurDisk;
    if (cosCurDisk < 0 && rmin > 0) rmincur /= cosHalf;
    rzmin[k].set(rmincur,rmin*sinCurDisk);

    G4double rmaxcur = rtor + rmax*cosCurDisk;
    if (cosCurDisk > 0) rmaxcur /= cosHalf;
    rzmax[k].set(rmaxcur,rmax*sinCurDisk);

    G4double sinTmpDisk = sinCurDisk;
    sinCurDisk = sinCurDisk*cosStepDisk + cosCurDisk*sinStepDisk;
    cosCurDisk = cosCurDisk*cosStepDisk - sinTmpDisk*sinStepDisk;
  }

  // Loop along slices in Phi. The extent is calculated as cumulative
  // extent of the slices
  pMin =  kInfinity;
  pMax = -kInfinity;
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);
  G4double sinCur1 = 0, cosCur1 = 0, sinCur2 = 0, cosCur2 = 0;
  for (G4int i=0; i<kphi+1; ++i)
  {
    if (i == 0)
    {
      sinCur1 = sinStart;
      cosCur1 = cosStart;
      sinCur2 = sinCur1*cosHalf + cosCur1*sinHalf;
      cosCur2 = cosCur1*cosHalf - sinCur1*sinHalf;
    }
    else
    {
      sinCur1 = sinCur2;
      cosCur1 = cosCur2;
      sinCur2 = (i == kphi) ? sinEnd : sinCur1*cosStep + cosCur1*sinStep;
      cosCur2 = (i == kphi) ? cosEnd : cosCur1*cosStep - sinCur1*sinStep;
    }
    for (G4int k=0; k<NDISK; ++k)
    {
      G4double r1 = rzmin[k].x(), r2 = rzmax[k].x();
      G4double z1 = rzmin[k].y(), z2 = rzmax[k].y();
      pols[k][0].set(r1*cosCur1,r1*sinCur1,z1);
      pols[k][1].set(r2*cosCur1,r2*sinCur1,z2);
      pols[k][2].set(r2*cosCur2,r2*sinCur2,z2);
      pols[k][3].set(r1*cosCur2,r1*sinCur2,z1);
    }
    pols[NDISK] = pols[0];

    // get bounding box of current slice
    G4TwoVector vmin,vmax;
    G4GeomTools::
      DiskExtent(rint,rext,sinCur1,cosCur1,sinCur2,cosCur2,vmin,vmax);
    bmin.setX(vmin.x()); bmin.setY(vmin.y());
    bmax.setX(vmax.x()); bmax.setY(vmax.y());

    // set bounding envelope for current slice and adjust extent
    G4double emin,emax;
    G4BoundingEnvelope benv(bmin,bmax,polygons);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax)) continue;
    if (emin < pMin) pMin = emin;
    if (emax > pMax) pMax = emax;
    if (eminlim > pMin && emaxlim < pMax) break; // max possible extent
  }
  return (pMin < pMax);
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization

G4Polyhedron* G4UTorus::CreatePolyhedron() const
{
  return new G4PolyhedronTorus(GetRmin(),
                               GetRmax(),
                               GetRtor(),
                               GetSPhi(),
                               GetDPhi());
}

#endif  // G4GEOM_USE_USOLIDS
