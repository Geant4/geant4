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
// Implementation for G4UPara wrapper class
//
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------

#include "G4Para.hh"
#include "G4UPara.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
//  Constructor - set & check half widths

G4UPara::G4UPara(const G4String& pName,
                       G4double pDx, G4double pDy, G4double pDz,
                       G4double pAlpha, G4double pTheta, G4double pPhi)
  : Base_t(pName, pDx, pDy, pDz, pAlpha, pTheta, pPhi)
{
  fTalpha = std::tan(pAlpha);
  fTthetaCphi = std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta)*std::sin(pPhi);
  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor - design of trapezoid based on 8 vertices

G4UPara::G4UPara( const G4String& pName,
                  const G4ThreeVector pt[8] )
  : Base_t(pName)
{
  // Find dimensions and trigonometric values
  //
  G4double fDx = (pt[3].x() - pt[2].x())*0.5;
  G4double fDy = (pt[2].y() - pt[1].y())*0.5;
  G4double fDz = pt[7].z();
  SetDimensions(fDx, fDy, fDz);
  CheckParameters(); // check dimensions

  fTalpha = (pt[2].x() + pt[3].x() - pt[1].x() - pt[0].x())*0.25/fDy;
  fTthetaCphi = (pt[4].x() + fDy*fTalpha + fDx)/fDz;
  fTthetaSphi = (pt[4].y() + fDy)/fDz;
  SetAlpha(std::atan(fTalpha));
  SetTheta(std::atan(std::sqrt(fTthetaSphi*fTthetaSphi
                              + fTthetaCphi*fTthetaCphi)));
  SetPhi (std::atan2(fTthetaSphi, fTthetaCphi));
  MakePlanes();

  // Recompute vertices
  //
  G4ThreeVector v[8];
  G4double DyTalpha = fDy*fTalpha;
  G4double DzTthetaSphi = fDz*fTthetaSphi;
  G4double DzTthetaCphi = fDz*fTthetaCphi;
  v[0].set(-DzTthetaCphi-DyTalpha-fDx, -DzTthetaSphi-fDy, -fDz);
  v[1].set(-DzTthetaCphi-DyTalpha+fDx, -DzTthetaSphi-fDy, -fDz);
  v[2].set(-DzTthetaCphi+DyTalpha-fDx, -DzTthetaSphi+fDy, -fDz);
  v[3].set(-DzTthetaCphi+DyTalpha+fDx, -DzTthetaSphi+fDy, -fDz);
  v[4].set( DzTthetaCphi-DyTalpha-fDx,  DzTthetaSphi-fDy,  fDz);
  v[5].set( DzTthetaCphi-DyTalpha+fDx,  DzTthetaSphi-fDy,  fDz);
  v[6].set( DzTthetaCphi+DyTalpha-fDx,  DzTthetaSphi+fDy,  fDz);
  v[7].set( DzTthetaCphi+DyTalpha+fDx,  DzTthetaSphi+fDy,  fDz);

  // Compare with original vertices
  //
  for (G4int i=0; i<8; ++i)
  {
    G4double delx = std::abs(pt[i].x() - v[i].x());
    G4double dely = std::abs(pt[i].y() - v[i].y());
    G4double delz = std::abs(pt[i].z() - v[i].z());
    G4double discrepancy = std::max(std::max(delx,dely),delz);
    if (discrepancy > 0.1*kCarTolerance)
    {
      std::ostringstream message;
      G4long oldprc = message.precision(16);
      message << "Invalid vertice coordinates for Solid: " << GetName()
              << "\nVertix #" << i << ", discrepancy = " << discrepancy
              << "\n  original   : " << pt[i]
              << "\n  recomputed : " << v[i];
      G4cout.precision(oldprc);
      G4Exception("G4UPara::G4UPara()", "GeomSolids0002",
                  FatalException, message);

    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency

G4UPara::G4UPara( __void__& a )
  : Base_t(a)
{
  SetAllParameters(1., 1., 1., 0., 0., 0.);
  fRebuildPolyhedron = false;
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UPara::~G4UPara()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UPara::G4UPara(const G4UPara& rhs)
  : Base_t(rhs), fTalpha(rhs.fTalpha),
    fTthetaCphi(rhs.fTthetaCphi),fTthetaSphi(rhs.fTthetaSphi)
{
  for (G4int i=0; i<4; ++i) { fPlanes[i] = rhs.fPlanes[i]; }
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UPara& G4UPara::operator = (const G4UPara& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   Base_t::operator=(rhs);

   // Copy data
   //
   fTalpha = rhs.fTalpha;
   fTthetaCphi = rhs.fTthetaCphi;
   fTthetaSphi = rhs.fTthetaSphi;
   for (G4int i=0; i<4; ++i) { fPlanes[i] = rhs.fPlanes[i]; }

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers

G4double G4UPara::GetZHalfLength() const
{
  return GetZ();
}
G4double G4UPara::GetYHalfLength() const
{
  return GetY();
}
G4double G4UPara::GetXHalfLength() const
{
  return GetX();
}
G4ThreeVector G4UPara::GetSymAxis() const
{
  return G4ThreeVector(fTthetaCphi,fTthetaSphi,1.).unit();
}
G4double G4UPara::GetTanAlpha() const
{
  return fTalpha;
}

G4double G4UPara::GetPhi() const       
{
   return std::atan2(fTthetaSphi,fTthetaCphi);
}

G4double G4UPara::GetTheta() const
{
   return std::atan(std::sqrt(fTthetaCphi*fTthetaCphi
                              +fTthetaSphi*fTthetaSphi));
}

G4double G4UPara::GetAlpha() const
{
  return std::atan(fTalpha);
}

void G4UPara::SetXHalfLength(G4double val)
{
  SetDimensions(val, GetY(), GetZ());
  fRebuildPolyhedron = true;

  CheckParameters();
  MakePlanes();
}
void G4UPara::SetYHalfLength(G4double val)
{
  SetDimensions(GetX(), val, GetZ());
  fRebuildPolyhedron = true;

  CheckParameters();
  MakePlanes();
}
void G4UPara::SetZHalfLength(G4double val)
{
  SetDimensions(GetX(), GetY(), val);
  fRebuildPolyhedron = true;

  CheckParameters();
  MakePlanes();
}
void G4UPara::SetAlpha(G4double alpha)
{
  Base_t::SetAlpha(alpha);
  fTalpha = std::tan(alpha);
  fRebuildPolyhedron = true;

  MakePlanes();
}
void G4UPara::SetTanAlpha(G4double val)
{
  fTalpha = val;
  fRebuildPolyhedron = true;

  MakePlanes();
}
void G4UPara::SetThetaAndPhi(double pTheta, double pPhi)
{
  Base_t::SetThetaAndPhi(pTheta, pPhi);
  G4double tanTheta = std::tan(pTheta);
  fTthetaCphi = tanTheta*std::cos(pPhi);
  fTthetaSphi = tanTheta*std::sin(pPhi);
  fRebuildPolyhedron = true;

  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - set and check half-widths

void G4UPara::SetAllParameters(G4double pDx, G4double pDy, G4double pDz,
                               G4double pAlpha, G4double pTheta, G4double pPhi)
{
  // Reset data of the base class
  fRebuildPolyhedron = true;

  // Set parameters
  SetDimensions(pDx, pDy, pDz);
  Base_t::SetAlpha(pAlpha);
  Base_t::SetThetaAndPhi(pTheta, pPhi);
  fTalpha = std::tan(pAlpha);
  fTthetaCphi = std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta)*std::sin(pPhi);

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Check dimensions

void G4UPara::CheckParameters()
{
  if (GetX() < 2*kCarTolerance ||
      GetY() < 2*kCarTolerance ||
      GetZ() < 2*kCarTolerance)
  {
    std::ostringstream message;
    message << "Invalid (too small or negative) dimensions for Solid: "
            << GetName()
            << "\n  X - " << GetX()
            << "\n  Y - " << GetY()
            << "\n  Z - " << GetZ();
    G4Exception("G4UPara::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Set side planes

void G4UPara::MakePlanes()
{
  G4ThreeVector vx(1, 0, 0);
  G4ThreeVector vy(fTalpha, 1, 0);
  G4ThreeVector vz(fTthetaCphi, fTthetaSphi, 1);

  // Set -Y & +Y planes
  //
  G4ThreeVector ynorm = (vx.cross(vz)).unit();

  fPlanes[0].a = 0.;
  fPlanes[0].b = ynorm.y();
  fPlanes[0].c = ynorm.z();
  fPlanes[0].d = fPlanes[0].b*GetY(); // point (0,fDy,0) is on plane

  fPlanes[1].a =  0.;
  fPlanes[1].b = -fPlanes[0].b;
  fPlanes[1].c = -fPlanes[0].c;
  fPlanes[1].d =  fPlanes[0].d;

  // Set -X & +X planes
  //
  G4ThreeVector xnorm = (vz.cross(vy)).unit();

  fPlanes[2].a = xnorm.x();
  fPlanes[2].b = xnorm.y();
  fPlanes[2].c = xnorm.z();
  fPlanes[2].d = fPlanes[2].a*GetZ(); // point (fDx,0,0) is on plane

  fPlanes[3].a = -fPlanes[2].a;
  fPlanes[3].b = -fPlanes[2].b;
  fPlanes[3].c = -fPlanes[2].c;
  fPlanes[3].d =  fPlanes[2].d;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification

void G4UPara::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*(G4Para*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UPara::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double dz = GetZHalfLength();
  G4double dx = GetXHalfLength();
  G4double dy = GetYHalfLength();

  G4double x0 = dz*fTthetaCphi;
  G4double x1 = dy*GetTanAlpha();
  G4double xmin =
    std::min(
    std::min(
    std::min(-x0-x1-dx,-x0+x1-dx),x0-x1-dx),x0+x1-dx);
  G4double xmax =
    std::max(
    std::max(
    std::max(-x0-x1+dx,-x0+x1+dx),x0-x1+dx),x0+x1+dx);

  G4double y0 = dz*fTthetaSphi;
  G4double ymin = std::min(-y0-dy,y0-dy);
  G4double ymax = std::max(-y0+dy,y0+dy);

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
    G4Exception("G4UPara::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4UPara::CalculateExtent( const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                       G4double& pMin, G4double& pMax ) const
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
  G4double dz = GetZHalfLength();
  G4double dx = GetXHalfLength();
  G4double dy = GetYHalfLength();

  G4double x0 = dz*fTthetaCphi;
  G4double x1 = dy*GetTanAlpha();
  G4double y0 = dz*fTthetaSphi;

  G4ThreeVectorList baseA(4), baseB(4);
  baseA[0].set(-x0-x1-dx,-y0-dy,-dz);
  baseA[1].set(-x0-x1+dx,-y0-dy,-dz);
  baseA[2].set(-x0+x1+dx,-y0+dy,-dz);
  baseA[3].set(-x0+x1-dx,-y0+dy,-dz);

  baseB[0].set(+x0-x1-dx, y0-dy, dz);
  baseB[1].set(+x0-x1+dx, y0-dy, dz);
  baseB[2].set(+x0+x1+dx, y0+dy, dz);
  baseB[3].set(+x0+x1-dx, y0+dy, dz);

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4UPara::Clone() const
{
  return new G4UPara(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

G4Polyhedron* G4UPara::CreatePolyhedron () const
{
  return new G4PolyhedronPara(GetX(), GetY(), GetZ(),
                              GetAlpha(), GetTheta(), GetPhi());
}

#endif  // G4GEOM_USE_USOLIDS
