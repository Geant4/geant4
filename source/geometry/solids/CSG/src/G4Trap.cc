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
// $Id: G4Trap.cc 104561 2017-06-06 07:54:54Z gcosmo $
//
// class G4Trap
//
// Implementation for G4Trap class
//
// History:
//
// 18.04.17 E.Tcherniaev: complete revision, speed-up
// 23.09.16 E.Tcherniaev: use G4BoundingEnvelope for CalculateExtent(),
//                      removed CreateRotatedVertices()
// 28.04.05 V.Grichine: new SurfaceNormal according to J. Apostolakis proposal 
// 26.04.05 V.Grichine: new SurfaceNormal is default 
// 19.04.05 V.Grichine: bug fixed in G4Trap("name",G4ThreeVector[8] vp)
// 12.12.04 V.Grichine: SurfaceNormal with edges/vertices 
// 15.11.04 V.Grichine: bug fixed in G4Trap("name",G4ThreeVector[8] vp)
// 13.12.99 V.Grichine: bug fixed in DistanceToIn(p,v)
// 19.11.99 V.Grichine: kUndef was added to Eside enum
// 04.06.99 S.Giani: Fixed CalculateExtent in rotated case. 
// 08.12.97 J.Allison: Added "nominal" constructor and method SetAllParameters.
// 01.11.96 V.Grichine: Costructor for Right Angular Wedge from STEP, G4Trd/Para
// 09.09.96 V.Grichine: Final modifications before to commit
// 21.03.95 P.Kent: Modified for `tolerant' geometry
//
///////////////////////////////////////////////////////////////////////////////

#include "G4Trap.hh"

#if !defined(G4GEOM_USE_UTRAP)

#include "globals.hh"
#include "G4GeomTools.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

using namespace CLHEP;

//////////////////////////////////////////////////////////////////////////
//
// Constructor - check and set half-widths as well as angles: 
// final check of coplanarity

G4Trap::G4Trap( const G4String& pName,
                      G4double pDz,
                      G4double pTheta, G4double pPhi,
                      G4double pDy1, G4double pDx1, G4double pDx2,
                      G4double pAlp1,
                      G4double pDy2, G4double pDx3, G4double pDx4,
                      G4double pAlp2)
  : G4CSGSolid(pName), halfCarTolerance(0.5*kCarTolerance)
{
  fDz = pDz;
  fTthetaCphi = std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta)*std::sin(pPhi);

  fDy1 = pDy1; fDx1 = pDx1; fDx2 = pDx2; fTalpha1 = std::tan(pAlp1);
  fDy2 = pDy2; fDx3 = pDx3; fDx4 = pDx4; fTalpha2 = std::tan(pAlp2);

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor - Design of trapezoid based on 8 G4ThreeVector parameters, 
// which are its vertices. Checking of planarity with preparation of 
// fPlanes[] and than calculation of other members

G4Trap::G4Trap( const G4String& pName,
                const G4ThreeVector pt[8] )
  : G4CSGSolid(pName), halfCarTolerance(0.5*kCarTolerance)
{
  // Start with check of centering - the center of gravity trap line
  // should cross the origin of frame
  //
  if (!(   pt[0].z() < 0
        && pt[0].z() == pt[1].z()
        && pt[0].z() == pt[2].z()
        && pt[0].z() == pt[3].z()

        && pt[4].z() > 0
        && pt[4].z() == pt[5].z()
        && pt[4].z() == pt[6].z()
        && pt[4].z() == pt[7].z()

        && std::fabs( pt[0].z() + pt[4].z() ) < kCarTolerance

        && pt[0].y() == pt[1].y()
        && pt[2].y() == pt[3].y()
        && pt[4].y() == pt[5].y()
        && pt[6].y() == pt[7].y()

        && std::fabs(pt[0].y()+pt[2].y()+pt[4].y()+pt[6].y()) < kCarTolerance
        && std::fabs(pt[0].x()+pt[1].x()+pt[4].x()+pt[5].x() +
                     pt[2].x()+pt[3].x()+pt[6].x()+pt[7].x()) < kCarTolerance ))
  {
    std::ostringstream message;
    message << "Invalid vertice coordinates for Solid: " << GetName();
    G4Exception("G4Trap::G4Trap()", "GeomSolids0002",
                FatalException, message);
  }
    
  // Set parameters
  //
  fDz = (pt[7]).z();
      
  fDy1     = ((pt[2]).y()-(pt[1]).y())*0.5;
  fDx1     = ((pt[1]).x()-(pt[0]).x())*0.5;
  fDx2     = ((pt[3]).x()-(pt[2]).x())*0.5;
  fTalpha1 = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/fDy1;

  fDy2     = ((pt[6]).y()-(pt[5]).y())*0.5;
  fDx3     = ((pt[5]).x()-(pt[4]).x())*0.5;
  fDx4     = ((pt[7]).x()-(pt[6]).x())*0.5;
  fTalpha2 = ((pt[6]).x()+(pt[7]).x()-(pt[5]).x()-(pt[4]).x())*0.25/fDy2;

  fTthetaCphi = ((pt[4]).x()+fDy2*fTalpha2+fDx3)/fDz;
  fTthetaSphi = ((pt[4]).y()+fDy2)/fDz;

  CheckParameters();
  MakePlanes(pt);
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor for Right Angular Wedge from STEP

G4Trap::G4Trap( const G4String& pName,
                      G4double pZ,
                      G4double pY,
                      G4double pX, G4double pLTX )
  : G4CSGSolid(pName), halfCarTolerance(0.5*kCarTolerance)
{
  fDz  = 0.5*pZ; fTthetaCphi = 0; fTthetaSphi = 0;
  fDy1 = 0.5*pY; fDx1 = 0.5*pX; fDx2 = 0.5*pLTX; fTalpha1 = 0.5*(pLTX - pX)/pY;
  fDy2 = fDy1;   fDx3 = fDx1;   fDx4 = fDx2;     fTalpha2 = fTalpha1;

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor for G4Trd

G4Trap::G4Trap( const G4String& pName,
                      G4double pDx1,  G4double pDx2,
                      G4double pDy1,  G4double pDy2,
                      G4double pDz )
  : G4CSGSolid(pName), halfCarTolerance(0.5*kCarTolerance), fTrapType(0)
{
  fDz  = pDz;  fTthetaCphi = 0; fTthetaSphi = 0;
  fDy1 = pDy1; fDx1 = pDx1; fDx2 = pDx1; fTalpha1 = 0;
  fDy2 = pDy2; fDx3 = pDx2; fDx4 = pDx2; fTalpha2 = 0;

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor for G4Para

G4Trap::G4Trap( const G4String& pName,
                      G4double pDx, G4double pDy,
                      G4double pDz,
                      G4double pAlpha,
                      G4double pTheta, G4double pPhi )
  : G4CSGSolid(pName), halfCarTolerance(0.5*kCarTolerance)
{
  fDz = pDz;
  fTthetaCphi = std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta)*std::sin(pPhi);

  fDy1 = pDy; fDx1 = pDx; fDx2 = pDx; fTalpha1 = std::tan(pAlpha);
  fDy2 = pDy; fDx3 = pDx; fDx4 = pDx; fTalpha2 = fTalpha1;

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Nominal constructor for G4Trap whose parameters are to be set by
// a G4VParamaterisation later.  Check and set half-widths as well as
// angles: final check of coplanarity

G4Trap::G4Trap( const G4String& pName )
  : G4CSGSolid (pName), halfCarTolerance(0.5*kCarTolerance),
    fDz(1.), fTthetaCphi(0.), fTthetaSphi(0.),
    fDy1(1.), fDx1(1.), fDx2(1.), fTalpha1(0.),
    fDy2(1.), fDx3(1.), fDx4(1.), fTalpha2(0.)
{
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Trap::G4Trap( __void__& a )
  : G4CSGSolid(a), halfCarTolerance(0.5*kCarTolerance),
    fDz(1.), fTthetaCphi(0.), fTthetaSphi(0.),
    fDy1(1.), fDx1(1.), fDx2(1.), fTalpha1(0.),
    fDy2(1.), fDx3(1.), fDx4(1.), fTalpha2(0.)
{
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Trap::~G4Trap()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Trap::G4Trap(const G4Trap& rhs)
  : G4CSGSolid(rhs), halfCarTolerance(rhs.halfCarTolerance),
    fDz(rhs.fDz), fTthetaCphi(rhs.fTthetaCphi), fTthetaSphi(rhs.fTthetaSphi),
    fDy1(rhs.fDy1), fDx1(rhs.fDx1), fDx2(rhs.fDx2), fTalpha1(rhs.fTalpha1),
    fDy2(rhs.fDy2), fDx3(rhs.fDx3), fDx4(rhs.fDx4), fTalpha2(rhs.fTalpha2)
{
  for (G4int i=0; i<4; ++i) { fPlanes[i] = rhs.fPlanes[i]; }
  fTrapType = rhs.fTrapType;
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Trap& G4Trap::operator = (const G4Trap& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4CSGSolid::operator=(rhs);

  // Copy data
  //
  halfCarTolerance = rhs.halfCarTolerance;
  fDz = rhs.fDz; fTthetaCphi = rhs.fTthetaCphi; fTthetaSphi = rhs.fTthetaSphi;
  fDy1 = rhs.fDy1; fDx1 = rhs.fDx1; fDx2 = rhs.fDx2; fTalpha1 = rhs.fTalpha1;
  fDy2 = rhs.fDy2; fDx3 = rhs.fDx3; fDx4 = rhs.fDx4; fTalpha2 = rhs.fTalpha2;
  for (G4int i=0; i<4; ++i) { fPlanes[i] = rhs.fPlanes[i]; }
  fTrapType = rhs.fTrapType;
  return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - check and set half-widths
// as well as angles: final check of coplanarity

void G4Trap::SetAllParameters ( G4double pDz,
                                G4double pTheta,
                                G4double pPhi,
                                G4double pDy1,
                                G4double pDx1,
                                G4double pDx2,
                                G4double pAlp1,
                                G4double pDy2,
                                G4double pDx3,
                                G4double pDx4,
                                G4double pAlp2 )
{
  // Reset data of the base class
  fCubicVolume = 0;
  fSurfaceArea = 0;
  fRebuildPolyhedron = true;

  // Set parameters
  fDz = pDz;
  fTthetaCphi = std::tan(pTheta)*std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta)*std::sin(pPhi);

  fDy1 = pDy1; fDx1 = pDx1; fDx2 = pDx2; fTalpha1 = std::tan(pAlp1);
  fDy2 = pDy2; fDx3 = pDx3; fDx4 = pDx4; fTalpha2 = std::tan(pAlp2);

  CheckParameters();
  MakePlanes();
}

//////////////////////////////////////////////////////////////////////////
//
// Check length parameters

void G4Trap::CheckParameters()
{
  if (fDz<=0 ||
      fDy1<=0 || fDx1<=0 || fDx2<=0 ||
      fDy2<=0 || fDx3<=0 || fDx4<=0)
  {
    std::ostringstream message;
    message << "Invalid Length Parameters for Solid: " << GetName()
            << "\n  X - " <<fDx1<<", "<<fDx2<<", "<<fDx3<<", "<<fDx4
            << "\n  Y - " <<fDy1<<", "<<fDy2
            << "\n  Z - " <<fDz;
    G4Exception("G4Trap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Compute vertices and set side planes

void G4Trap::MakePlanes()
{
  G4double DzTthetaCphi = fDz*fTthetaCphi;
  G4double DzTthetaSphi = fDz*fTthetaSphi;
  G4double Dy1Talpha1   = fDy1*fTalpha1;
  G4double Dy2Talpha2   = fDy2*fTalpha2;

  G4ThreeVector pt[8] =
  {
    G4ThreeVector(-DzTthetaCphi-Dy1Talpha1-fDx1,-DzTthetaSphi-fDy1,-fDz),
    G4ThreeVector(-DzTthetaCphi-Dy1Talpha1+fDx1,-DzTthetaSphi-fDy1,-fDz),
    G4ThreeVector(-DzTthetaCphi+Dy1Talpha1-fDx2,-DzTthetaSphi+fDy1,-fDz),
    G4ThreeVector(-DzTthetaCphi+Dy1Talpha1+fDx2,-DzTthetaSphi+fDy1,-fDz),
    G4ThreeVector( DzTthetaCphi-Dy2Talpha2-fDx3, DzTthetaSphi-fDy2, fDz),
    G4ThreeVector( DzTthetaCphi-Dy2Talpha2+fDx3, DzTthetaSphi-fDy2, fDz),
    G4ThreeVector( DzTthetaCphi+Dy2Talpha2-fDx4, DzTthetaSphi+fDy2, fDz),
    G4ThreeVector( DzTthetaCphi+Dy2Talpha2+fDx4, DzTthetaSphi+fDy2, fDz)
  };

  MakePlanes(pt);
}

//////////////////////////////////////////////////////////////////////////
//
// Set side planes, check planarity

void G4Trap::MakePlanes(const G4ThreeVector pt[8])
{
  G4int iface[4][4] = { {0,4,5,1}, {2,3,7,6}, {0,2,6,4}, {1,5,7,3} };
  G4String side[4] = { "~-Y", "~+Y", "~-X", "~+X" };

  for (G4int i=0; i<4; ++i)
  {
    if (MakePlane(pt[iface[i][0]],
                  pt[iface[i][1]],
                  pt[iface[i][2]],
                  pt[iface[i][3]],
                  fPlanes[i])) continue;

    // Non planar side face
    G4ThreeVector normal(fPlanes[i].a,fPlanes[i].b,fPlanes[i].c);
    G4double dmax = 0;
    for (G4int k=0; k<4; ++k)
    {
      G4double dist = normal.dot(pt[iface[i][k]]) + fPlanes[i].d;
      if (std::abs(dist) > std::abs(dmax)) dmax = dist;
    }
    std::ostringstream message;
    message << "Side face " << side[i] << " is not planar for solid: "
            << GetName() << "\nDiscrepancy: " << dmax/mm << " mm\n";
    StreamInfo(message);
    G4Exception("G4Trap::MakePlanes()", "GeomSolids0002",
                FatalException, message);
  }

  // Define type of trapezoid
  fTrapType = 0;
  if (fPlanes[0].b  == -1 && fPlanes[1].b == 1 &&
      std::abs(fPlanes[0].a) < DBL_EPSILON &&
      std::abs(fPlanes[0].c) < DBL_EPSILON &&
      std::abs(fPlanes[1].a) < DBL_EPSILON &&
      std::abs(fPlanes[1].c) < DBL_EPSILON)
  {
    fTrapType = 1; // YZ section is a rectangle ...
    if (std::abs(fPlanes[2].a + fPlanes[3].a) < DBL_EPSILON &&
        std::abs(fPlanes[2].c - fPlanes[3].c) < DBL_EPSILON &&
        fPlanes[2].b == 0 &&
        fPlanes[3].b == 0)
    {
      fTrapType = 2; // ... and XZ section is a isosceles trapezoid
      fPlanes[2].a = -fPlanes[3].a;
      fPlanes[2].c =  fPlanes[3].c;
    }
  }
}

///////////////////////////////////////////////////////////////////////
//
// Calculate the coef's of the plane p1->p2->p3->p4->p1
// where the ThreeVectors 1-4 are in anti-clockwise order when viewed
// from infront of the plane (i.e. from normal direction).
//
// Return true if the points are coplanar, false otherwise

G4bool G4Trap::MakePlane( const G4ThreeVector& p1,
                          const G4ThreeVector& p2,
                          const G4ThreeVector& p3,
                          const G4ThreeVector& p4,
                                TrapSidePlane& plane )
{
  G4ThreeVector normal = ((p4 - p2).cross(p3 - p1)).unit();
  if (std::abs(normal.x()) < DBL_EPSILON) normal.setX(0); 
  if (std::abs(normal.y()) < DBL_EPSILON) normal.setY(0); 
  if (std::abs(normal.z()) < DBL_EPSILON) normal.setZ(0); 
  normal = normal.unit();

  G4ThreeVector centre = (p1 + p2 + p3 + p4)*0.25;
  plane.a =  normal.x();
  plane.b =  normal.y();
  plane.c =  normal.z();
  plane.d = -normal.dot(centre);

  // compute distances and check planarity
  G4double d1 = std::abs(normal.dot(p1) + plane.d);
  G4double d2 = std::abs(normal.dot(p2) + plane.d);
  G4double d3 = std::abs(normal.dot(p3) + plane.d);
  G4double d4 = std::abs(normal.dot(p4) + plane.d);
  G4double dmax = std::max(std::max(std::max(d1,d2),d3),d4);
  
  return (dmax > 1000 * kCarTolerance) ? false : true;
}

///////////////////////////////////////////////////////////////////////
//
// Get volume

G4double G4Trap::GetCubicVolume()
{
  if (fCubicVolume == 0)
  {
    G4ThreeVector pt[8];
    GetVertices(pt);
 
    G4double dz  = pt[4].z() - pt[0].z();
    G4double dy1 = pt[2].y() - pt[0].y();
    G4double dx1 = pt[1].x() - pt[0].x();
    G4double dx2 = pt[3].x() - pt[2].x();
    G4double dy2 = pt[6].y() - pt[4].y();
    G4double dx3 = pt[5].x() - pt[4].x();
    G4double dx4 = pt[7].x() - pt[6].x();

    fCubicVolume = ((dx1 + dx2 + dx3 + dx4)*(dy1 + dy2) +
                    (dx4 + dx3 - dx2 - dx1)*(dy2 - dy1)/3)*dz*0.125;
  }
  return fCubicVolume;
}

///////////////////////////////////////////////////////////////////////
//
// Get surface area

G4double G4Trap::GetSurfaceArea()
{
  if (fSurfaceArea == 0)
  {
    G4ThreeVector pt[8];
    G4int iface [6][4] =
      { {0,1,3,2}, {0,4,5,1}, {2,3,7,6}, {0,2,6,4}, {1,5,7,3}, {4,6,7,5} };

    GetVertices(pt);
    for (G4int i=0; i<6; ++i)
    {
      fSurfaceArea += G4GeomTools::QuadAreaNormal(pt[iface[i][0]],
                                                  pt[iface[i][1]],
                                                  pt[iface[i][2]],
                                                  pt[iface[i][3]]).mag();
    }
  }
  return fSurfaceArea;
}

///////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Trap::ComputeDimensions(       G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

///////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Trap::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4ThreeVector pt[8];
  GetVertices(pt);

  G4double xmin = kInfinity, xmax = -kInfinity;
  G4double ymin = kInfinity, ymax = -kInfinity;
  for (G4int i=0; i<8; ++i)
  {
    G4double x = pt[i].x();
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    G4double y = pt[i].y();
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
  }

  G4double dz   = GetZHalfLength();
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
    G4Exception("G4Trap::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

///////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Trap::CalculateExtent( const EAxis pAxis,
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
  G4ThreeVector pt[8];
  GetVertices(pt);

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

///////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Trap::Inside( const G4ThreeVector& p ) const
{
  if (fTrapType == 2)      // YZ section is a rectangle and
  {                        // XZ section is an isosceles trapezoid
    G4double dy = std::max(std::abs(p.z())-fDz,std::abs(p.y())+fPlanes[1].d);
    G4double dx = fPlanes[3].a*std::abs(p.x())+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,dx);

    if (dist > halfCarTolerance) return kOutside;
    return (dist > -halfCarTolerance) ? kSurface : kInside;
  }
  else if (fTrapType == 1) // YZ section is a rectangle
  {
    G4double dy = std::max(std::abs(p.z())-fDz,std::abs(p.y())+fPlanes[1].d);
    G4double dx1 = fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    G4double dx2 = fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,std::max(dx1,dx2));

    if (dist > halfCarTolerance) return kOutside;
    return (dist > -halfCarTolerance) ? kSurface : kInside;
  }
  else                     // General case
  {
    G4double dz = std::abs(p.z())-fDz;
    G4double dy1 = fPlanes[0].b*p.y()+fPlanes[0].c*p.z()+fPlanes[0].d;
    G4double dy2 = fPlanes[1].b*p.y()+fPlanes[1].c*p.z()+fPlanes[1].d;
    G4double dy = std::max(dz,std::max(dy1,dy2));

    G4double dx1 = fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    G4double dx2 = fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,std::max(dx1,dx2));

    if (dist > halfCarTolerance) return kOutside;
    return (dist > -halfCarTolerance) ? kSurface : kInside;
  }
}

///////////////////////////////////////////////////////////////////////
//
// Determine side, and return corresponding normal

G4ThreeVector G4Trap::SurfaceNormal( const G4ThreeVector& p ) const
{
  // Check Z faces
  //
  G4double nz = 0;
  if (std::abs(std::abs(p.z()) - fDz) <= halfCarTolerance)
  {
    nz = (p.z() < 0) ? -1 : 1;
  }

  // Check Y faces
  //
  G4double ny = 0;
  if (fTrapType > 0)   // YZ section is a rectangle
  {
    G4double dist = std::abs(p.y()) + fPlanes[1].d;
    if (std::abs(dist) <= halfCarTolerance) ny = (p.y() < 0) ? -1 : 1;
  }
  else
  {
    for (G4int i=0; i<2; ++i)
    {
      G4double dist = fPlanes[i].b*p.y() + fPlanes[i].c*p.z() + fPlanes[i].d;
      if (std::abs(dist) > halfCarTolerance) continue;
      ny  = fPlanes[i].b;
      nz += fPlanes[i].c;
      break;
    }
  }

  // Check X faces
  //
  G4double nx = 0;
  if (fTrapType == 2)   // YZ section is a rectangle and
  {                     // XZ section is an isosceles trapezoid
    G4double dist = fPlanes[3].a*std::abs(p.x())
                  + fPlanes[3].c*p.z() + fPlanes[3].d;
    if (std::abs(dist) <= halfCarTolerance)
    {
      nx  = (p.x() < 0) ? -fPlanes[3].a : fPlanes[3].a;
      nz += fPlanes[3].c;
    }
  }
  else
  {
    for (G4int i=2; i<4; ++i)
    {
      G4double dist = fPlanes[i].a*p.x() +
                      fPlanes[i].b*p.y() + fPlanes[i].c*p.z() + fPlanes[i].d;
      if (std::abs(dist) > halfCarTolerance) continue;
      nx  = fPlanes[i].a;
      ny += fPlanes[i].b;
      nz += fPlanes[i].c;
      break;
    }
  }

  // Return normal
  //
  G4int nsurf = nx*nx + ny*ny + nz*nz + 0.5;                  // get magnitude
  if (nsurf == 1)      return G4ThreeVector(nx,ny,nz);
  else if (nsurf != 0) return G4ThreeVector(nx,ny,nz).unit(); // edge or corner
  else
  {
    // Point is not on the surface
    //
#ifdef G4CSGDEBUG
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is not on surface (!?) of solid: "
            << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc) ;
    G4Exception("G4Trap::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Trap::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  G4double dist = -DBL_MAX;
  G4int iside = 0;
  for (G4int i=0; i<4; ++i)
  {
    G4double d = fPlanes[i].a*p.x() +
                 fPlanes[i].b*p.y() +
                 fPlanes[i].c*p.z() + fPlanes[i].d;
    if (d > dist) { dist = d; iside = i; }
  }

  G4double distz = std::abs(p.z()) - fDz;
  if (dist > distz)
    return G4ThreeVector(fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c);
  else
    return G4ThreeVector(0, 0, (p.z() < 0) ? -1 : 1);
}

///////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside
//  - return kInfinity if no intersection

G4double G4Trap::DistanceToIn(const G4ThreeVector& p,
                              const G4ThreeVector& v ) const
{
  // Z intersections
  //
  if ((std::abs(p.z()) - fDz) >= -halfCarTolerance && p.z()*v.z() >= 0)
    return kInfinity;
  G4double invz = (-v.z() == 0) ? DBL_MAX : -1./v.z();
  G4double dz = (invz < 0) ? fDz : -fDz; 
  G4double tzmin = (p.z() + dz)*invz;
  G4double tzmax = (p.z() - dz)*invz;

  // Y intersections
  //
  G4double tymin = 0, tymax = DBL_MAX;
  G4int i = 0;
  for ( ; i<2; ++i)
  { 
    G4double cosa = fPlanes[i].b*v.y() + fPlanes[i].c*v.z();
    G4double dist = fPlanes[i].b*p.y() + fPlanes[i].c*p.z() + fPlanes[i].d;
    if (dist >= -halfCarTolerance)
    {
      if (cosa >= 0) return kInfinity;
      G4double tmp  = -dist/cosa;
      if (tymin < tmp) tymin = tmp;
    }
    else if (cosa > 0)
    {
      G4double tmp  = -dist/cosa;
      if (tymax > tmp) tymax = tmp;
    } 
  }

  // Z intersections
  //
  G4double txmin = 0, txmax = DBL_MAX;
  for ( ; i<4; ++i)
  { 
    G4double cosa = fPlanes[i].a*v.x()+fPlanes[i].b*v.y()+fPlanes[i].c*v.z();
    G4double dist = fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z() +
                    fPlanes[i].d;
    if (dist >= -halfCarTolerance)
    {
      if (cosa >= 0) return kInfinity;
      G4double tmp  = -dist/cosa;
      if (txmin < tmp) txmin = tmp;
    }
    else if (cosa > 0)
    {
      G4double tmp  = -dist/cosa;
      if (txmax > tmp) txmax = tmp;
    } 
  }

  // Find distance
  //
  G4double tmin = std::max(std::max(txmin,tymin),tzmin);
  G4double tmax = std::min(std::min(txmax,tymax),tzmax);
     
  if (tmax <= tmin + halfCarTolerance) return kInfinity; // touch or no hit
  return (tmin < halfCarTolerance ) ? 0. : tmin;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// This is the best fast estimation of the shortest distance to trap
// - Returns 0 is ThreeVector inside

G4double G4Trap::DistanceToIn( const G4ThreeVector& p ) const
{
  if (fTrapType == 2)      // YZ section is a rectangle and
  {                        // XZ section is an isosceles trapezoid
    G4double dy = std::max(std::abs(p.z())-fDz,std::abs(p.y())+fPlanes[1].d);
    G4double dx = fPlanes[3].a*std::abs(p.x())+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,dx);

    return (dist > 0) ? dist : 0.;
  }
  else if (fTrapType == 1) // YZ section is a rectangle
  {
    G4double dy = std::max(std::abs(p.z())-fDz,std::abs(p.y())+fPlanes[1].d);
    G4double dx1 = fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    G4double dx2 = fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,std::max(dx1,dx2));

    return (dist > 0) ? dist : 0.;
  }
  else                     // General case
  {
    G4double dz = std::abs(p.z())-fDz;
    G4double dy1 = fPlanes[0].b*p.y()+fPlanes[0].c*p.z()+fPlanes[0].d;
    G4double dy2 = fPlanes[1].b*p.y()+fPlanes[1].c*p.z()+fPlanes[1].d;
    G4double dy = std::max(dz,std::max(dy1,dy2));

    G4double dx1 = fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    G4double dx2 = fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,std::max(dx1,dx2));

    return (dist > 0) ? dist : 0.;
  }
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside and
// find normal at exit point, if required
// - when leaving the surface, return 0

G4double G4Trap::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool *validNorm, G4ThreeVector *n) const
{
  // Z intersections
  //
  if ((std::abs(p.z()) - fDz) >= -halfCarTolerance && p.z()*v.z() > 0)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0, 0, (p.z() < 0) ? -1 : 1);
    }
    return 0;
  }
  G4double vz = v.z();
  G4double tmax = (vz == 0) ? DBL_MAX : (std::copysign(fDz,vz) - p.z())/vz;
  G4int iside = (vz < 0) ? -4 : -2; // little trick: (-4+3)=-1, (-2+3)=+1

  // Y intersections
  //
  G4int i = 0;
  for ( ; i<2; ++i)
  {
    G4double cosa = fPlanes[i].b*v.y() + fPlanes[i].c*v.z();
    if (cosa > 0)
    {
      G4double dist = fPlanes[i].b*p.y() + fPlanes[i].c*p.z() + fPlanes[i].d;
      if (dist >= -halfCarTolerance)
      {
        if (calcNorm)
        {
          *validNorm = true;
          n->set(0, fPlanes[i].b, fPlanes[i].c);
        }
        return 0;
      }
      G4double tmp = -dist/cosa;
      if (tmax > tmp) { tmax = tmp; iside = i; }
    }
  }

  // X intersections
  //
  for ( ; i<4; ++i)
  {
    G4double cosa = fPlanes[i].a*v.x()+fPlanes[i].b*v.y()+fPlanes[i].c*v.z();
    if (cosa > 0)
    {
      G4double dist = fPlanes[i].a*p.x() +
                      fPlanes[i].b*p.y() + fPlanes[i].c*p.z() + fPlanes[i].d;
      if (dist >= -halfCarTolerance)
      {
        if (calcNorm)
        {
           *validNorm = true;
           n->set(fPlanes[i].a, fPlanes[i].b, fPlanes[i].c);
        }
        return 0;
      }
      G4double tmp = -dist/cosa;
      if (tmax > tmp) { tmax = tmp; iside = i; }
    }
  }

  // Set normal, if required, and return distance
  //
  if (calcNorm) 
  {
    *validNorm = true;
    if (iside < 0)
      n->set(0, 0, iside + 3); // (-4+3)=-1, (-2+3)=+1
    else
      n->set(fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c);
  }
  return tmax;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is ThreeVector outside

G4double G4Trap::DistanceToOut( const G4ThreeVector& p ) const
{
#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Point p is outside (!?) of solid: " << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc) ;
    G4Exception("G4Trap::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
  }
#endif
  if (fTrapType == 2)      // YZ section is a rectangle and
  {                        // XZ section is an isosceles trapezoid
    G4double dy = std::max(std::abs(p.z())-fDz,std::abs(p.y())+fPlanes[1].d);
    G4double dx = fPlanes[3].a*std::abs(p.x())+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,dx);

    return (dist < 0) ? -dist : 0.;
  }
  else if (fTrapType == 1) // YZ section is a rectangle
  {
    G4double dy = std::max(std::abs(p.z())-fDz,std::abs(p.y())+fPlanes[1].d);
    G4double dx1 = fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    G4double dx2 = fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,std::max(dx1,dx2));

    return (dist < 0) ? -dist : 0.;
  }
  else                     // General case
  {
    G4double dz = std::abs(p.z())-fDz;
    G4double dy1 = fPlanes[0].b*p.y()+fPlanes[0].c*p.z()+fPlanes[0].d;
    G4double dy2 = fPlanes[1].b*p.y()+fPlanes[1].c*p.z()+fPlanes[1].d;
    G4double dy = std::max(dz,std::max(dy1,dy2));

    G4double dx1 = fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    G4double dx2 = fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    G4double dist = std::max(dy,std::max(dx1,dx2));

    return (dist < 0) ? -dist : 0.;
  }
}

////////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Trap::GetEntityType() const
{
  return G4String("G4Trap");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Trap::Clone() const
{
  return new G4Trap(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Trap::StreamInfo( std::ostream& os ) const
{
  G4double phi    = std::atan2(fTthetaSphi,fTthetaCphi);
  G4double theta  = std::atan(std::sqrt(fTthetaCphi*fTthetaCphi
                                       +fTthetaSphi*fTthetaSphi));
  G4double alpha1 = std::atan(fTalpha1);
  G4double alpha2 = std::atan(fTalpha2);
  G4String signDegree = "\u00B0"; 

  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid: " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Trap\n"
     << " Parameters:\n"
     << "    half length Z: " << fDz/mm << " mm\n"
     << "    half length Y, face -Dz: " << fDy1/mm << " mm\n"
     << "    half length X, face -Dz, side -Dy1: " << fDx1/mm << " mm\n"
     << "    half length X, face -Dz, side +Dy1: " << fDx2/mm << " mm\n"
     << "    half length Y, face +Dz: " << fDy2/mm << " mm\n"
     << "    half length X, face +Dz, side -Dy2: " << fDx3/mm << " mm\n"
     << "    half length X, face +Dz, side +Dy2: " << fDx4/mm << " mm\n"
     << "    theta: " << theta/degree << signDegree << "\n"
     << "    phi: " << phi/degree << signDegree << "\n"
     << "    alpha, face -Dz: " << alpha1/degree << signDegree << "\n"
     << "    alpha, face +Dz: " << alpha2/degree << signDegree << "\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Compute vertices from planes

void G4Trap::GetVertices(G4ThreeVector pt[8]) const
{
  for (G4int i=0; i<8; ++i)
  {
    G4int iy = (i==0 || i==1 || i==4 || i==5) ? 0 : 1;
    G4int ix = (i==0 || i==2 || i==4 || i==6) ? 2 : 3;
    G4double z = (i < 4) ? -fDz : fDz;
    G4double y = -(fPlanes[iy].c*z + fPlanes[iy].d)/fPlanes[iy].b;
    G4double x = -(fPlanes[ix].b*y + fPlanes[ix].c*z
                   + fPlanes[ix].d)/fPlanes[ix].a;
    pt[i].set(x,y,z);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface

G4ThreeVector G4Trap::GetPointOnSurface() const
{
  G4ThreeVector pt[8];
  G4int iface [6][4] =
    { {0,1,3,2}, {0,4,5,1}, {2,3,7,6}, {0,2,6,4}, {1,5,7,3}, {4,6,7,5} }; 
  G4double sface[6];

  GetVertices(pt);
  G4double stotal = 0;
  for (G4int i=0; i<6; ++i)
  {
    G4double ss = G4GeomTools::QuadAreaNormal(pt[iface[i][0]],
                                              pt[iface[i][1]],
                                              pt[iface[i][2]],
                                              pt[iface[i][3]]).mag();
    stotal  += ss;
    sface[i] = stotal;
  }

  // Select face
  //
  G4double select = stotal*G4UniformRand();
  G4int k = 5;
  if (select <= sface[4]) k = 4;
  if (select <= sface[3]) k = 3;
  if (select <= sface[2]) k = 2;
  if (select <= sface[1]) k = 1;
  if (select <= sface[0]) k = 0;

  // Select sub-triangle
  //
  G4int i0 = iface[k][0];
  G4int i1 = iface[k][1];
  G4int i2 = iface[k][2];
  G4int i3 = iface[k][3];
  G4double s1 = G4GeomTools::TriangleAreaNormal(pt[i0],pt[i1],pt[i3]).mag();
  G4double s2 = G4GeomTools::TriangleAreaNormal(pt[i2],pt[i1],pt[i3]).mag();
  if ((s1+s2)*G4UniformRand() > s1) i0 = i2;

  // Generate point
  //
  G4double u = G4UniformRand();
  G4double v = G4UniformRand();
  if (u + v > 1.) { u = 1. - u; v = 1. - v; }
  return (1.-u-v)*pt[i0] + u*pt[i1] + v*pt[i3];
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Trap::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Trap::CreatePolyhedron () const
{
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);
  G4double alpha1 = std::atan(fTalpha1);
  G4double alpha2 = std::atan(fTalpha2);
  G4double theta = std::atan(std::sqrt(fTthetaCphi*fTthetaCphi
                                      +fTthetaSphi*fTthetaSphi));

  return new G4PolyhedronTrap(fDz, theta, phi,
                              fDy1, fDx1, fDx2, alpha1,
                              fDy2, fDx3, fDx4, alpha2);
}

#endif
