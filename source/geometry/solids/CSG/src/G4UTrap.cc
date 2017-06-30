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
// Implementation for G4UTrap wrapper class
// --------------------------------------------------------------------

#include "G4Trap.hh"
#include "G4UTrap.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructors
//
G4UTrap::G4UTrap( const G4String& pName,
                        G4double pdz,
                        G4double pTheta, G4double pPhi,
                        G4double pdy1, G4double pdx1, G4double pdx2,
                        G4double pAlp1,
                        G4double pdy2, G4double pdx3, G4double pdx4,
                        G4double pAlp2 )
  : Base_t(pName, pdz, pTheta, pPhi, pdy1, pdx1, pdx2,
           pAlp1, pdy2, pdx3, pdx4, pAlp2)
{
}

G4UTrap::G4UTrap( const G4String& pName,
                  const G4ThreeVector pt[8] )
  : Base_t(pName)
{
  SetPlanes(pt);
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pZ,
                        G4double pY,
                        G4double pX, G4double pLTX )
  : Base_t(pName, pZ, pY, pX, pLTX)
{
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pdx1,  G4double pdx2,
                        G4double pdy1,  G4double pdy2,
                        G4double pdz )
  : Base_t(pName, pdx1, pdx2, pdy1, pdy2, pdz)
{
}

G4UTrap::G4UTrap(const G4String& pName,
                       G4double pdx, G4double pdy, G4double pdz,
                       G4double pAlpha, G4double pTheta, G4double pPhi )
  : Base_t(pName, pdx, pdy, pdz, pAlpha, pTheta, pPhi)
{
}

G4UTrap::G4UTrap( const G4String& pName )
  : Base_t(pName)
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTrap::G4UTrap( __void__& a )
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTrap::~G4UTrap()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTrap::G4UTrap(const G4UTrap& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTrap& G4UTrap::operator = (const G4UTrap& rhs) 
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

G4double G4UTrap::GetZHalfLength() const
{
  return GetDz();
}
G4double G4UTrap::GetYHalfLength1() const
{
  return GetDy1();
}
G4double G4UTrap::GetXHalfLength1() const
{
  return GetDx1();
}
G4double G4UTrap::GetXHalfLength2() const
{
  return GetDx2();
}
G4double G4UTrap::GetYHalfLength2() const
{
  return GetDy2();
}
G4double G4UTrap::GetXHalfLength3() const
{
  return GetDx3();
}
G4double G4UTrap::GetXHalfLength4() const
{
  return GetDx4();
}
G4double G4UTrap::GetThetaCphi() const
{
  return GetTanThetaCosPhi();
}
G4double G4UTrap::GetThetaSphi() const
{
  return GetTanThetaSinPhi();
}
TrapSidePlane G4UTrap::GetSidePlane(G4int n) const
{
  TrapSidePlane plane;
  plane.a = GetStruct().GetPlane(n).fA;
  plane.b = GetStruct().GetPlane(n).fB;
  plane.c = GetStruct().GetPlane(n).fC;
  plane.d = GetStruct().GetPlane(n).fD;
  return plane;
}
G4ThreeVector G4UTrap::GetSymAxis() const
{
  G4double tanThetaSphi = GetTanThetaSinPhi();
  G4double tanThetaCphi = GetTanThetaCosPhi();
  G4double tan2Theta = tanThetaSphi*tanThetaSphi + tanThetaCphi*tanThetaCphi;
  G4double cosTheta = 1.0 / std::sqrt(1 + tan2Theta);
  return G4ThreeVector(tanThetaCphi*cosTheta, tanThetaSphi*cosTheta, cosTheta);
}

void G4UTrap::SetAllParameters(G4double pDz, G4double pTheta, G4double pPhi,
                               G4double pDy1, G4double pDx1, G4double pDx2,
                               G4double pAlp1,
                               G4double pDy2, G4double pDx3, G4double pDx4,
                               G4double pAlp2)
{
  SetDz(pDz);
  SetDy1(pDy1);
  SetDy2(pDy2);
  SetDx1(pDx1);
  SetDx2(pDx2);
  SetDx3(pDx3);
  SetDx4(pDx4);
  SetTanAlpha1(std::tan(pAlp1));
  SetTanAlpha1(std::tan(pAlp2));
  // last two will also reset cached variables
  SetTheta(pTheta);
  SetPhi(pPhi);
  fRebuildPolyhedron = true;
}

void G4UTrap::SetPlanes(const G4ThreeVector pt[8])
{
  U3Vector upt[8];
  for (unsigned int i=0; i<8; ++i)
  {
    upt[i] = U3Vector(pt[i].x(), pt[i].y(), pt[i].z());
  }
  fromCornersToParameters(upt);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4UTrap::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Trap*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4UTrap::Clone() const
{
  return new G4UTrap(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTrap::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  TrapSidePlane planes[4];
  for (G4int i=0; i<4; ++i) { planes[i] = GetSidePlane(i); }

  G4double xmin = kInfinity, xmax = -kInfinity;
  G4double ymin = kInfinity, ymax = -kInfinity;
  G4double dz   = GetZHalfLength();
  for (G4int i=0; i<8; ++i)
  {
    G4int iy = (i==0 || i==1 || i==4 || i==5) ? 0 : 1;
    G4int ix = (i==0 || i==2 || i==4 || i==6) ? 2 : 3;
    G4double z = (i < 4) ? -dz : dz;
    G4double y = -(planes[iy].c*z + planes[iy].d)/planes[iy].b;
    G4double x = -(planes[ix].b*y + planes[ix].c*z + planes[ix].d)/planes[ix].a;
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }

  pMin.set(xmin,ymin,-dz);
  pMax.set(xmax,ymax, dz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UTrap::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  // Check consistency of bounding boxes
  //
  if (checkBBox)
  {
    G4double tolerance = kCarTolerance;
    U3Vector vmin, vmax;
    Extent(vmin,vmax);
    if (std::abs(pMin.x()-vmin.x()) > tolerance ||
        std::abs(pMin.y()-vmin.y()) > tolerance ||
        std::abs(pMin.z()-vmin.z()) > tolerance ||
        std::abs(pMax.x()-vmax.x()) > tolerance ||
        std::abs(pMax.y()-vmax.y()) > tolerance ||
        std::abs(pMax.z()-vmax.z()) > tolerance)
    {
      std::ostringstream message;
      message << "Inconsistency in bounding boxes for solid: "
              << GetName() << " !"
              << "\nBBox min: wrapper = " << pMin << " solid = " << vmin
              << "\nBBox max: wrapper = " << pMax << " solid = " << vmax;
      G4Exception("G4UTrap::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UTrap::CalculateExtent(const EAxis pAxis,
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
  TrapSidePlane planes[4];
  for (G4int i=0; i<4; ++i) { planes[i] = GetSidePlane(i); }

  G4ThreeVector pt[8];
  G4double dz = GetZHalfLength();
  for (G4int i=0; i<8; ++i)
  {
    G4int iy = (i==0 || i==1 || i==4 || i==5) ? 0 : 1;
    G4int ix = (i==0 || i==2 || i==4 || i==6) ? 2 : 3;
    G4double z = (i < 4) ? -dz : dz;
    G4double y = -(planes[iy].c*z + planes[iy].d)/planes[iy].b;
    G4double x = -(planes[ix].b*y + planes[ix].c*z + planes[ix].d)/planes[ix].a;
    pt[i].set(x,y,z);
  }

  G4ThreeVectorList baseA(4), baseB(4);
  baseA[0] = pt[0];
  baseA[1] = pt[1];
  baseA[2] = pt[3];
  baseA[3] = pt[2];

  baseB[0] = pt[4];
  baseB[1] = pt[5];
  baseB[2] = pt[7];
  baseB[3] = pt[6];

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
G4Polyhedron* G4UTrap::CreatePolyhedron() const
{
  G4double fTthetaSphi = GetThetaSphi();
  G4double fTthetaCphi = GetThetaCphi();
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);
  G4double alpha1 = std::atan(GetTanAlpha1());
  G4double alpha2 = std::atan(GetTanAlpha2());
  G4double theta = std::atan(std::sqrt(fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi));

  return new G4PolyhedronTrap(GetZHalfLength(), theta, phi,
                              GetYHalfLength1(),
                              GetXHalfLength1(), GetXHalfLength2(), alpha1,
                              GetYHalfLength2(),
                              GetXHalfLength3(), GetXHalfLength4(), alpha2);
}

#endif  // G4GEOM_USE_USOLIDS
