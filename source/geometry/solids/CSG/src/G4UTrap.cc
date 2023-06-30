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
// Implementation for G4UTrap wrapper class
//
// 13.09.13 G.Cosmo, CERN/PH
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
  G4ThreeVector pt[8];
  CheckParameters();
  GetVertices(pt);
  CheckPlanarity(pt);
}

G4UTrap::G4UTrap( const G4String& pName,
                  const G4ThreeVector pt[8] )
  : Base_t(pName)
{
  // Start with check of centering - the center of gravity trap line
  // should cross the origin of frame
  if (  pt[0].z() >= 0
        || pt[0].z() != pt[1].z()
        || pt[0].z() != pt[2].z()
        || pt[0].z() != pt[3].z()

        || pt[4].z() <= 0
        || pt[4].z() != pt[5].z()
        || pt[4].z() != pt[6].z()
        || pt[4].z() != pt[7].z()

        || std::abs( pt[0].z() + pt[4].z() ) >= kCarTolerance

        || pt[0].y() != pt[1].y()
        || pt[2].y() != pt[3].y()
        || pt[4].y() != pt[5].y()
        || pt[6].y() != pt[7].y()

        || std::abs(pt[0].y()+pt[2].y()+pt[4].y()+pt[6].y()) >= kCarTolerance
        || std::abs(pt[0].x()+pt[1].x()+pt[4].x()+pt[5].x() +
                    pt[2].x()+pt[3].x()+pt[6].x()+pt[7].x()) >= kCarTolerance )
  {
    std::ostringstream message;
    message << "Invalid vertice coordinates for Solid: " << GetName();
    G4Exception("G4UTrap::G4UTrap()", "GeomSolids0002",
                FatalException, message);
  }

  SetPlanes(pt);
  CheckParameters();
  CheckPlanarity(pt);
}

// Constructor for Right Angular Wedge from STEP; this is different from
// the UnplacedTrapezoid constructor taking four double's and constructing
// a Trd1.
G4UTrap::G4UTrap( const G4String& pName,
                        G4double pZ,
                        G4double pY,
                        G4double pX, G4double pLTX )
  : Base_t(pName, 0.5*pZ, /*theta=*/0, /*phi=*/0, /*dy1=*/0.5*pY,
           /*dx1=*/0.5*pX, /*dx2=*/0.5*pLTX, /*Alpha1=*/0.5*(pLTX - pX)/pY,
           /*dy2=*/0.5*pY, /*dx3=*/0.5*pX, /*dx4=*/0.5*pLTX,
           /*Alpha2=*/0.5*(pLTX - pX)/pY)
{
  CheckParameters();
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pdx1,  G4double pdx2,
                        G4double pdy1,  G4double pdy2,
                        G4double pdz )
  : Base_t(pName, pdx1, pdx2, pdy1, pdy2, pdz)
{
  CheckParameters();
}

G4UTrap::G4UTrap(const G4String& pName,
                       G4double pdx, G4double pdy, G4double pdz,
                       G4double pAlpha, G4double pTheta, G4double pPhi )
  : Base_t(pName, pdx, pdy, pdz, pAlpha, pTheta, pPhi)
{
  CheckParameters();
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
G4UTrap::~G4UTrap() = default;

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
// Accessors
//
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
G4double G4UTrap::GetTanAlpha1() const
{
  return Base_t::GetTanAlpha1();
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
G4double G4UTrap::GetTanAlpha2() const
{
  return Base_t::GetTanAlpha2();
}
G4double G4UTrap::GetPhi() const       
{
  return Base_t::GetPhi();
}
G4double G4UTrap::GetTheta() const
{
  return Base_t::GetTheta();
}
G4double G4UTrap::GetAlpha1() const
{
  return Base_t::GetAlpha1();
}
G4double G4UTrap::GetAlpha2() const
{
  return Base_t::GetAlpha2();
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
  return {tanThetaCphi*cosTheta, tanThetaSphi*cosTheta, cosTheta};
}

//////////////////////////////////////////////////////////////////////////
//
// Modifier
//
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

  G4ThreeVector pt[8];
  CheckParameters();
  GetVertices(pt);
  CheckPlanarity(pt);
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters using eight vertices
//
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
// Check dimensions
//
void G4UTrap::CheckParameters() const
{
  G4double fDz  = GetZHalfLength();
  G4double fDy1 = GetYHalfLength1();
  G4double fDx1 = GetXHalfLength1();
  G4double fDx2 = GetXHalfLength2();
  G4double fDy2 = GetYHalfLength2();
  G4double fDx3 = GetXHalfLength3();
  G4double fDx4 = GetXHalfLength4();

  if (fDz<=0 ||
      fDy1<=0 || fDx1<=0 || fDx2<=0 ||
      fDy2<=0 || fDx3<=0 || fDx4<=0)
  {
    std::ostringstream message;
    message << "Invalid Length Parameters for Solid: " << GetName()
            << "\n  X - " <<fDx1<<", "<<fDx2<<", "<<fDx3<<", "<<fDx4
            << "\n  Y - " <<fDy1<<", "<<fDy2
            << "\n  Z - " <<fDz;
    G4Exception("G4UTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Compute coordinates of vertices
//
void G4UTrap::GetVertices(G4ThreeVector pt[8]) const
{
  G4double fDz      = GetZHalfLength();
  G4double fDy1     = GetYHalfLength1();
  G4double fDx1     = GetXHalfLength1();
  G4double fDx2     = GetXHalfLength2();
  G4double fDy2     = GetYHalfLength2();
  G4double fDx3     = GetXHalfLength3();
  G4double fDx4     = GetXHalfLength4();
  G4double phi      = GetPhi();
  G4double theta    = GetTheta();
  G4double fTalpha1 = GetTanAlpha1();
  G4double fTalpha2 = GetTanAlpha2();

  G4double DzTthetaCphi = fDz*std::tan(theta)*std::cos(phi);
  G4double DzTthetaSphi = fDz*std::tan(theta)*std::sin(phi);
  G4double Dy1Talpha1   = fDy1*fTalpha1;
  G4double Dy2Talpha2   = fDy2*fTalpha2;

  pt[0].set(-DzTthetaCphi-Dy1Talpha1-fDx1,-DzTthetaSphi-fDy1,-fDz);
  pt[1].set(-DzTthetaCphi-Dy1Talpha1+fDx1,-DzTthetaSphi-fDy1,-fDz);
  pt[2].set(-DzTthetaCphi+Dy1Talpha1-fDx2,-DzTthetaSphi+fDy1,-fDz);
  pt[3].set(-DzTthetaCphi+Dy1Talpha1+fDx2,-DzTthetaSphi+fDy1,-fDz);
  pt[4].set( DzTthetaCphi-Dy2Talpha2-fDx3, DzTthetaSphi-fDy2, fDz);
  pt[5].set( DzTthetaCphi-Dy2Talpha2+fDx3, DzTthetaSphi-fDy2, fDz);
  pt[6].set( DzTthetaCphi+Dy2Talpha2-fDx4, DzTthetaSphi+fDy2, fDz);
  pt[7].set( DzTthetaCphi+Dy2Talpha2+fDx4, DzTthetaSphi+fDy2, fDz);
}

/////////////////////////////////////////////////////////////////////////
//
// Check planarity of lateral planes
//
void G4UTrap::CheckPlanarity(const G4ThreeVector pt[8]) const
{
  constexpr G4int iface[4][4] = { {0,4,5,1}, {2,3,7,6}, {0,2,6,4}, {1,5,7,3} };
  const static G4String side[4] = { "~-Y", "~+Y", "~-X", "~+X" };

  for (G4int i=0; i<4; ++i)
  {
    TrapSidePlane plane = GetSidePlane(i);
    G4double dmax = 0;
    for (G4int k=0; k<4; ++k)
    {
      const G4ThreeVector p = pt[iface[i][k]];
      G4double dist = plane.a*p.x() + plane.b*p.y() + plane.c*p.z() + plane.d;
      if (std::abs(dist) > std::abs(dmax)) dmax = dist;
    }
    if (std::abs(dmax) > 1000 * kCarTolerance)
    {
      std::ostringstream message;
      message << "Side face " << side[i] << " is not planar for solid: "
              << GetName() << "\nDiscrepancy: " << dmax/mm << " mm\n";
      StreamInfo(message);
      G4Exception("G4UTrap::CheckPlanarity()", "GeomSolids0002",
                  FatalException, message);
    }
  }
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
    return exist = pMin < pMax;
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
  return new G4PolyhedronTrap(GetZHalfLength(), GetTheta(), GetPhi(),
                              GetYHalfLength1(),
                              GetXHalfLength1(), GetXHalfLength2(), GetAlpha1(),
                              GetYHalfLength2(),
                              GetXHalfLength3(), GetXHalfLength4(), GetAlpha2());
}

#endif  // G4GEOM_USE_USOLIDS
