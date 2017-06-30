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
// $Id: G4Cons.cc 104316 2017-05-24 13:04:23Z gcosmo $
// GEANT4 tag $Name: $
//
//
// class G4Cons
//
// Implementation for G4Cons class
//
// History:
//
// 03.10.16 E.Tcherniaev: use G4BoundingEnvelope for CalculateExtent(),
//                      removed CreateRotatedVertices()
// 04.09.14 T.Nikitina: Fix typo error in GetPointOnSurface() when 
//                      GetRadiusInRing() was introduced
//                      Fix DistanceToIn(p,v) for points on the Surface,
//                      error was reported by OpticalEscape test
// 05.04.12 M.Kelsey:   GetPointOnSurface() throw flat in sqrt(r)
// 12.10.09 T.Nikitina: Added to DistanceToIn(p,v) check on the direction in
//                      case of point on surface
// 03.05.05 V.Grichine: SurfaceNormal(p) according to J. Apostolakis proposal
// 13.09.96 V.Grichine: Review and final modifications
// ~1994    P.Kent: Created, as main part of the geometry prototype
// --------------------------------------------------------------------

#include "G4Cons.hh"

#if !defined(G4GEOM_USE_UCONS)

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"

using namespace CLHEP;
 
////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

// used by normal

enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

//////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

G4Cons::G4Cons( const G4String& pName,
                      G4double  pRmin1, G4double pRmax1,
                      G4double  pRmin2, G4double pRmax2,
                      G4double pDz,
                      G4double pSPhi, G4double pDPhi)
  : G4CSGSolid(pName), fRmin1(pRmin1), fRmin2(pRmin2),
    fRmax1(pRmax1), fRmax2(pRmax2), fDz(pDz), fSPhi(0.), fDPhi(0.)
{
  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  halfCarTolerance=kCarTolerance*0.5;
  halfRadTolerance=kRadTolerance*0.5;
  halfAngTolerance=kAngTolerance*0.5;

  // Check z-len
  //
  if ( pDz < 0 )
  {
    std::ostringstream message;
    message << "Invalid Z half-length for Solid: " << GetName() << G4endl
            << "        hZ = " << pDz;
    G4Exception("G4Cons::G4Cons()", "GeomSolids0002",
                FatalException, message);
  }

  // Check radii
  //
  if (((pRmin1>=pRmax1) || (pRmin2>=pRmax2) || (pRmin1<0)) && (pRmin2<0))
  {
    std::ostringstream message;
    message << "Invalid values of radii for Solid: " << GetName() << G4endl
            << "        pRmin1 = " << pRmin1 << ", pRmin2 = " << pRmin2
            << ", pRmax1 = " << pRmax1 << ", pRmax2 = " << pRmax2;
    G4Exception("G4Cons::G4Cons()", "GeomSolids0002",
                FatalException, message) ;
  }
  if( (pRmin1 == 0.0) && (pRmin2 > 0.0) ) { fRmin1 = 1e3*kRadTolerance ; }
  if( (pRmin2 == 0.0) && (pRmin1 > 0.0) ) { fRmin2 = 1e3*kRadTolerance ; }

  // Check angles
  //
  CheckPhiAngles(pSPhi, pDPhi);
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Cons::G4Cons( __void__& a )
  : G4CSGSolid(a), kRadTolerance(0.), kAngTolerance(0.),
    fRmin1(0.), fRmin2(0.), fRmax1(0.), fRmax2(0.), fDz(0.),
    fSPhi(0.), fDPhi(0.), sinCPhi(0.), cosCPhi(0.), cosHDPhiOT(0.),
    cosHDPhiIT(0.), sinSPhi(0.), cosSPhi(0.), sinEPhi(0.), cosEPhi(0.),
    fPhiFullCone(false), halfCarTolerance(0.), halfRadTolerance(0.),
    halfAngTolerance(0.)
{
}

///////////////////////////////////////////////////////////////////////
//
// Destructor

G4Cons::~G4Cons()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Cons::G4Cons(const G4Cons& rhs)
  : G4CSGSolid(rhs), kRadTolerance(rhs.kRadTolerance),
    kAngTolerance(rhs.kAngTolerance), fRmin1(rhs.fRmin1), fRmin2(rhs.fRmin2),
    fRmax1(rhs.fRmax1), fRmax2(rhs.fRmax2), fDz(rhs.fDz), fSPhi(rhs.fSPhi),
    fDPhi(rhs.fDPhi), sinCPhi(rhs.sinCPhi), cosCPhi(rhs.cosCPhi),
    cosHDPhiOT(rhs.cosHDPhiOT), cosHDPhiIT(rhs.cosHDPhiIT),
    sinSPhi(rhs.sinSPhi), cosSPhi(rhs.cosSPhi), sinEPhi(rhs.sinEPhi),
    cosEPhi(rhs.cosEPhi), fPhiFullCone(rhs.fPhiFullCone),
    halfCarTolerance(rhs.halfCarTolerance),
    halfRadTolerance(rhs.halfRadTolerance),
    halfAngTolerance(rhs.halfAngTolerance)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Cons& G4Cons::operator = (const G4Cons& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   kRadTolerance = rhs.kRadTolerance;
   kAngTolerance = rhs.kAngTolerance;
   fRmin1 = rhs.fRmin1; fRmin2 = rhs.fRmin2;
   fRmax1 = rhs.fRmax1; fRmax2 = rhs.fRmax2;
   fDz = rhs.fDz; fSPhi = rhs.fSPhi; fDPhi = rhs.fDPhi;
   sinCPhi = rhs.sinCPhi; cosCPhi = rhs.cosCPhi;
   cosHDPhiOT = rhs.cosHDPhiOT; cosHDPhiIT = rhs.cosHDPhiIT;
   sinSPhi = rhs.sinSPhi; cosSPhi = rhs.cosSPhi;
   sinEPhi = rhs.sinEPhi; cosEPhi = rhs.cosEPhi;
   fPhiFullCone = rhs.fPhiFullCone;
   halfCarTolerance = rhs.halfCarTolerance;
   halfRadTolerance = rhs.halfRadTolerance;
   halfAngTolerance = rhs.halfAngTolerance;

   return *this;
}

/////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Cons::Inside(const G4ThreeVector& p) const
{
  G4double r2, rl, rh, pPhi, tolRMin, tolRMax; // rh2, rl2 ;
  EInside in;

  if (std::fabs(p.z()) > fDz + halfCarTolerance )  { return in = kOutside; }
  else if(std::fabs(p.z()) >= fDz - halfCarTolerance )    { in = kSurface; }
  else                                                    { in = kInside;  }

  r2 = p.x()*p.x() + p.y()*p.y() ;
  rl = 0.5*(fRmin2*(p.z() + fDz) + fRmin1*(fDz - p.z()))/fDz ;
  rh = 0.5*(fRmax2*(p.z()+fDz)+fRmax1*(fDz-p.z()))/fDz;

  // rh2 = rh*rh;

  tolRMin = rl - halfRadTolerance;
  if ( tolRMin < 0 )  { tolRMin = 0; }
  tolRMax = rh + halfRadTolerance;

  if ( (r2<tolRMin*tolRMin) || (r2>tolRMax*tolRMax) ) { return in = kOutside; }

  if (rl) { tolRMin = rl + halfRadTolerance; }
  else    { tolRMin = 0.0; }
  tolRMax = rh - halfRadTolerance;
      
  if (in == kInside) // else it's kSurface already
  {
     if ( (r2 < tolRMin*tolRMin) || (r2 >= tolRMax*tolRMax) ) { in = kSurface; }
  }
  if ( !fPhiFullCone && ((p.x() != 0.0) || (p.y() != 0.0)) )
  {
    pPhi = std::atan2(p.y(),p.x()) ;

    if ( pPhi < fSPhi - halfAngTolerance  )             { pPhi += twopi; }
    else if ( pPhi > fSPhi + fDPhi + halfAngTolerance ) { pPhi -= twopi; }
    
    if ( (pPhi < fSPhi - halfAngTolerance) ||          
         (pPhi > fSPhi + fDPhi + halfAngTolerance) )  { return in = kOutside; }
      
    else if (in == kInside)  // else it's kSurface anyway already
    {
       if ( (pPhi < fSPhi + halfAngTolerance) || 
            (pPhi > fSPhi + fDPhi - halfAngTolerance) )  { in = kSurface; }
    }
  }
  else if ( !fPhiFullCone )  { in = kSurface; }

  return in ;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Cons::ComputeDimensions(      G4VPVParameterisation* p,
                               const G4int                  n,
                               const G4VPhysicalVolume*     pRep    )
{
  p->ComputeDimensions(*this,n,pRep) ;
}

///////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Cons::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
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
    G4Exception("G4Cons::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

///////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Cons::CalculateExtent( const EAxis              pAxis,
                                const G4VoxelLimits&     pVoxelLimit,
                                const G4AffineTransform& pTransform,
                                      G4double&          pMin,
                                      G4double&          pMax ) const
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

////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Cons::SurfaceNormal( const G4ThreeVector& p) const
{
  G4int noSurfaces = 0;
  G4double rho, pPhi;
  G4double distZ, distRMin, distRMax;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;
  G4double tanRMin, secRMin, pRMin, widRMin;
  G4double tanRMax, secRMax, pRMax, widRMax;

  G4ThreeVector norm, sumnorm(0.,0.,0.), nZ = G4ThreeVector(0.,0.,1.);
  G4ThreeVector nR, nr(0.,0.,0.), nPs, nPe;

  distZ = std::fabs(std::fabs(p.z()) - fDz);
  rho   = std::sqrt(p.x()*p.x() + p.y()*p.y());

  tanRMin  = (fRmin2 - fRmin1)*0.5/fDz;
  secRMin  = std::sqrt(1 + tanRMin*tanRMin);
  pRMin    = rho - p.z()*tanRMin;
  widRMin  = fRmin2 - fDz*tanRMin;
  distRMin = std::fabs(pRMin - widRMin)/secRMin;

  tanRMax  = (fRmax2 - fRmax1)*0.5/fDz;
  secRMax  = std::sqrt(1+tanRMax*tanRMax);
  pRMax    = rho - p.z()*tanRMax;
  widRMax  = fRmax2 - fDz*tanRMax;
  distRMax = std::fabs(pRMax - widRMax)/secRMax;

  if (!fPhiFullCone)   // Protected against (0,0,z) 
  {
    if ( rho )
    {
      pPhi = std::atan2(p.y(),p.x());

      if (pPhi  < fSPhi-halfCarTolerance)            { pPhi += twopi; }
      else if (pPhi > fSPhi+fDPhi+halfCarTolerance)  { pPhi -= twopi; }

      distSPhi = std::fabs( pPhi - fSPhi ); 
      distEPhi = std::fabs( pPhi - fSPhi - fDPhi ); 
    }
    else if( !(fRmin1) || !(fRmin2) )
    {
      distSPhi = 0.; 
      distEPhi = 0.; 
    }
    nPs = G4ThreeVector(std::sin(fSPhi), -std::cos(fSPhi), 0);
    nPe = G4ThreeVector(-std::sin(fSPhi+fDPhi), std::cos(fSPhi+fDPhi), 0);
  }
  if ( rho > halfCarTolerance )   
  {
    nR = G4ThreeVector(p.x()/rho/secRMax, p.y()/rho/secRMax, -tanRMax/secRMax);
    if (fRmin1 || fRmin2)
    {
      nr = G4ThreeVector(-p.x()/rho/secRMin,-p.y()/rho/secRMin,tanRMin/secRMin);
    }
  }

  if( distRMax <= halfCarTolerance )
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  if( (fRmin1 || fRmin2) && (distRMin <= halfCarTolerance) )
  {
    noSurfaces ++;
    sumnorm += nr;
  }
  if( !fPhiFullCone )   
  {
    if (distSPhi <= halfAngTolerance)
    {
      noSurfaces ++;
      sumnorm += nPs;
    }
    if (distEPhi <= halfAngTolerance) 
    {
      noSurfaces ++;
      sumnorm += nPe;
    }
  }
  if (distZ <= halfCarTolerance)  
  {
    noSurfaces ++;
    if ( p.z() >= 0.)  { sumnorm += nZ; }
    else               { sumnorm -= nZ; }
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4Cons::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" );
#endif 
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = sumnorm.unit(); }

  return norm ;
}

////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Cons::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  ENorm side ;
  G4ThreeVector norm ;
  G4double rho, phi ;
  G4double distZ, distRMin, distRMax, distSPhi, distEPhi, distMin ;
  G4double tanRMin, secRMin, pRMin, widRMin ;
  G4double tanRMax, secRMax, pRMax, widRMax ;

  distZ = std::fabs(std::fabs(p.z()) - fDz) ;
  rho   = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

  tanRMin  = (fRmin2 - fRmin1)*0.5/fDz ;
  secRMin  = std::sqrt(1 + tanRMin*tanRMin) ;
  pRMin    = rho - p.z()*tanRMin ;
  widRMin  = fRmin2 - fDz*tanRMin ;
  distRMin = std::fabs(pRMin - widRMin)/secRMin ;

  tanRMax  = (fRmax2 - fRmax1)*0.5/fDz ;
  secRMax  = std::sqrt(1+tanRMax*tanRMax) ;
  pRMax    = rho - p.z()*tanRMax ;
  widRMax  = fRmax2 - fDz*tanRMax ;
  distRMax = std::fabs(pRMax - widRMax)/secRMax ;
  
  if (distRMin < distRMax)  // First minimum
  {
    if (distZ < distRMin)
    {
      distMin = distZ ;
      side    = kNZ ;
    }
    else
    {
      distMin = distRMin ;
      side    = kNRMin ;
    }
  }
  else
  {
    if (distZ < distRMax)
    {
      distMin = distZ ;
      side    = kNZ ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax ;
    }
  }
  if ( !fPhiFullCone && rho )  // Protected against (0,0,z) 
  {
    phi = std::atan2(p.y(),p.x()) ;

    if (phi < 0)  { phi += twopi; }

    if (fSPhi < 0)  { distSPhi = std::fabs(phi - (fSPhi + twopi))*rho; }
    else            { distSPhi = std::fabs(phi - fSPhi)*rho; }

    distEPhi = std::fabs(phi - fSPhi - fDPhi)*rho ;

    // Find new minimum

    if (distSPhi < distEPhi)
    {
      if (distSPhi < distMin)  { side = kNSPhi; }
    }
    else 
    {
      if (distEPhi < distMin)  { side = kNEPhi; }
    }
  }    
  switch (side)
  {
    case kNRMin:      // Inner radius
      rho *= secRMin ;
      norm = G4ThreeVector(-p.x()/rho, -p.y()/rho, tanRMin/secRMin) ;
      break ;
    case kNRMax:      // Outer radius
      rho *= secRMax ;
      norm = G4ThreeVector(p.x()/rho, p.y()/rho, -tanRMax/secRMax) ;
      break ;
    case kNZ:         // +/- dz
      if (p.z() > 0)  { norm = G4ThreeVector(0,0,1);  }
      else            { norm = G4ThreeVector(0,0,-1); }
      break ;
    case kNSPhi:
      norm = G4ThreeVector(std::sin(fSPhi), -std::cos(fSPhi), 0) ;
      break ;
    case kNEPhi:
      norm=G4ThreeVector(-std::sin(fSPhi+fDPhi), std::cos(fSPhi+fDPhi), 0) ;
      break ;
    default:          // Should never reach this case...
      DumpInfo();
      G4Exception("G4Cons::ApproxSurfaceNormal()",
                  "GeomSolids1002", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break ;    
  }
  return norm ;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outside cone, compute intersection with rmax1*0.5
//        - if at valid phi,z return
//        - if inside outer cone, handle case when on tolerant outer cone
//          boundary and heading inwards(->0 to in)
//
// -> Compute intersection with inner cone, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - `if valid' implies tolerant checking of intersection points
// - z, phi intersection from Tubs

G4double G4Cons::DistanceToIn( const G4ThreeVector& p,
                               const G4ThreeVector& v   ) const
{
  G4double snxt = kInfinity ;      // snxt = default return value
  const G4double dRmax = 50*(fRmax1+fRmax2);// 100*(Rmax1+Rmax2)/2.

  G4double tanRMax,secRMax,rMaxAv,rMaxOAv ;  // Data for cones
  G4double tanRMin,secRMin,rMinAv,rMinOAv ;
  G4double rout,rin ;

  G4double tolORMin,tolORMin2,tolIRMin,tolIRMin2 ; // `generous' radii squared
  G4double tolORMax2,tolIRMax,tolIRMax2 ;
  G4double tolODz,tolIDz ;

  G4double Dist,sd,xi,yi,zi,ri=0.,risec,rhoi2,cosPsi ; // Intersection point vars

  G4double t1,t2,t3,b,c,d ;    // Quadratic solver variables 
  G4double nt1,nt2,nt3 ;
  G4double Comp ;

  G4ThreeVector Normal;

  // Cone Precalcs

  tanRMin = (fRmin2 - fRmin1)*0.5/fDz ;
  secRMin = std::sqrt(1.0 + tanRMin*tanRMin) ;
  rMinAv  = (fRmin1 + fRmin2)*0.5 ;

  if (rMinAv > halfRadTolerance)
  {
    rMinOAv = rMinAv - halfRadTolerance ;
  }
  else
  {
    rMinOAv = 0.0 ;
  }  
  tanRMax = (fRmax2 - fRmax1)*0.5/fDz ;
  secRMax = std::sqrt(1.0 + tanRMax*tanRMax) ;
  rMaxAv  = (fRmax1 + fRmax2)*0.5 ;
  rMaxOAv = rMaxAv + halfRadTolerance ;
   
  // Intersection with z-surfaces

  tolIDz = fDz - halfCarTolerance ;
  tolODz = fDz + halfCarTolerance ;

  if (std::fabs(p.z()) >= tolIDz)
  {
    if ( p.z()*v.z() < 0 )    // at +Z going in -Z or visa versa
    {
      sd = (std::fabs(p.z()) - fDz)/std::fabs(v.z()) ; // Z intersect distance

      if( sd < 0.0 )  { sd = 0.0; }                    // negative dist -> zero

      xi   = p.x() + sd*v.x() ;  // Intersection coords
      yi   = p.y() + sd*v.y() ;
      rhoi2 = xi*xi + yi*yi  ;

      // Check validity of intersection
      // Calculate (outer) tolerant radi^2 at intersecion

      if (v.z() > 0)
      {
        tolORMin  = fRmin1 - halfRadTolerance*secRMin ;
        tolIRMin  = fRmin1 + halfRadTolerance*secRMin ;
        tolIRMax  = fRmax1 - halfRadTolerance*secRMin ;
        // tolORMax2 = (fRmax1 + halfRadTolerance*secRMax)*
        //             (fRmax1 + halfRadTolerance*secRMax) ;
      }
      else
      {
        tolORMin  = fRmin2 - halfRadTolerance*secRMin ;
        tolIRMin  = fRmin2 + halfRadTolerance*secRMin ;
        tolIRMax  = fRmax2 - halfRadTolerance*secRMin ;
        // tolORMax2 = (fRmax2 + halfRadTolerance*secRMax)*
        //             (fRmax2 + halfRadTolerance*secRMax) ;
      }
      if ( tolORMin > 0 ) 
      {
        // tolORMin2 = tolORMin*tolORMin ;
        tolIRMin2 = tolIRMin*tolIRMin ;
      }
      else                
      {
        // tolORMin2 = 0.0 ;
        tolIRMin2 = 0.0 ;
      }
      if ( tolIRMax > 0 )  { tolIRMax2 = tolIRMax*tolIRMax; }     
      else                 { tolIRMax2 = 0.0; }
      
      if ( (tolIRMin2 <= rhoi2) && (rhoi2 <= tolIRMax2) )
      {
        if ( !fPhiFullCone && rhoi2 )
        {
          // Psi = angle made with central (average) phi of shape

          cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;

          if (cosPsi >= cosHDPhiIT)  { return sd; }
        }
        else
        {
          return sd;
        }
      }
    }
    else  // On/outside extent, and heading away  -> cannot intersect
    {
      return snxt ;  
    }
  }
    
// ----> Can not intersect z surfaces


// Intersection with outer cone (possible return) and
//                   inner cone (must also check phi)
//
// Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
//
// Intersects with x^2+y^2=(a*z+b)^2
//
// where a=tanRMax or tanRMin
//       b=rMaxAv  or rMinAv
//
// (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0 ;
//     t1                        t2                      t3  
//
//  \--------u-------/       \-----------v----------/ \---------w--------/
//

  t1   = 1.0 - v.z()*v.z() ;
  t2   = p.x()*v.x() + p.y()*v.y() ;
  t3   = p.x()*p.x() + p.y()*p.y() ;
  rin  = tanRMin*p.z() + rMinAv ;
  rout = tanRMax*p.z() + rMaxAv ;

  // Outer Cone Intersection
  // Must be outside/on outer cone for valid intersection

  nt1 = t1 - (tanRMax*v.z())*(tanRMax*v.z()) ;
  nt2 = t2 - tanRMax*v.z()*rout ;
  nt3 = t3 - rout*rout ;

  if (std::fabs(nt1) > kRadTolerance)  // Equation quadratic => 2 roots
  {
    b = nt2/nt1;
    c = nt3/nt1;
    d = b*b-c  ;
    if ( (nt3 > rout*rout*kRadTolerance*kRadTolerance*secRMax*secRMax)
      || (rout < 0) )
    {
      // If outside real cone (should be rho-rout>kRadTolerance*0.5
      // NOT rho^2 etc) saves a std::sqrt() at expense of accuracy

      if (d >= 0)
      {
          
        if ((rout < 0) && (nt3 <= 0))
        {
          // Inside `shadow cone' with -ve radius
          // -> 2nd root could be on real cone

          if (b>0) { sd = c/(-b-std::sqrt(d)); }
          else     { sd = -b + std::sqrt(d);   }
        }
        else
        {
          if ((b <= 0) && (c >= 0)) // both >=0, try smaller root
          {
            sd=c/(-b+std::sqrt(d));
          }
          else
          {
            if ( c <= 0 ) // second >=0
            {
              sd = -b + std::sqrt(d) ;
              if((sd<0) & (sd>-halfRadTolerance)) sd=0;
            }
            else  // both negative, travel away
            {
              return kInfinity ;
            }
          }
        }
        if ( sd >= 0 )  // If 'forwards'. Check z intersection
        {
          if ( sd>dRmax ) // Avoid rounding errors due to precision issues on
          {               // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd-std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
          zi = p.z() + sd*v.z() ;

          if (std::fabs(zi) <= tolODz)
          {
            // Z ok. Check phi intersection if reqd

            if ( fPhiFullCone )  { return sd; }
            else
            {
              xi     = p.x() + sd*v.x() ;
              yi     = p.y() + sd*v.y() ;
              ri     = rMaxAv + zi*tanRMax ;
              cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri ;

              if ( cosPsi >= cosHDPhiIT )  { return sd; }
            }
          }
        }                // end if (sd>0)
      }
    }
    else
    {
      // Inside outer cone
      // check not inside, and heading through G4Cons (-> 0 to in)

      if ( ( t3  > (rin + halfRadTolerance*secRMin)*
                   (rin + halfRadTolerance*secRMin) )
        && (nt2 < 0) && (d >= 0) && (std::fabs(p.z()) <= tolIDz) )
      {
        // Inside cones, delta r -ve, inside z extent
        // Point is on the Surface => check Direction using  Normal.dot(v)

        xi     = p.x() ;
        yi     = p.y()  ;
        risec  = std::sqrt(xi*xi + yi*yi)*secRMax ;
        Normal = G4ThreeVector(xi/risec,yi/risec,-tanRMax/secRMax) ;
        if ( !fPhiFullCone )
        {
          cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/std::sqrt(t3) ;
          if ( cosPsi >= cosHDPhiIT )
          {
            if ( Normal.dot(v) <= 0 )  { return 0.0; }
          }
        }
        else
        {             
          if ( Normal.dot(v) <= 0 )  { return 0.0; }
        }
      }
    }
  }
  else  //  Single root case 
  {
    if ( std::fabs(nt2) > kRadTolerance )
    {
      sd = -0.5*nt3/nt2 ;

      if ( sd < 0 )  { return kInfinity; }   // travel away
      else  // sd >= 0,  If 'forwards'. Check z intersection
      {
        zi = p.z() + sd*v.z() ;

        if ((std::fabs(zi) <= tolODz) && (nt2 < 0))
        {
          // Z ok. Check phi intersection if reqd

          if ( fPhiFullCone )  { return sd; }
          else
          {
            xi     = p.x() + sd*v.x() ;
            yi     = p.y() + sd*v.y() ;
            ri     = rMaxAv + zi*tanRMax ;
            cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri ;

            if (cosPsi >= cosHDPhiIT)  { return sd; }
          }
        }
      }
    }
    else  //    travel || cone surface from its origin
    {
      sd = kInfinity ;
    }
  }

  // Inner Cone Intersection
  // o Space is divided into 3 areas:
  //   1) Radius greater than real inner cone & imaginary cone & outside
  //      tolerance
  //   2) Radius less than inner or imaginary cone & outside tolarance
  //   3) Within tolerance of real or imaginary cones
  //      - Extra checks needed for 3's intersections
  //        => lots of duplicated code

  if (rMinAv)
  { 
    nt1 = t1 - (tanRMin*v.z())*(tanRMin*v.z()) ;
    nt2 = t2 - tanRMin*v.z()*rin ;
    nt3 = t3 - rin*rin ;
 
    if ( nt1 )
    {
      if ( nt3 > rin*kRadTolerance*secRMin )
      {
        // At radius greater than real & imaginary cones
        // -> 2nd root, with zi check

        b = nt2/nt1 ;
        c = nt3/nt1 ;
        d = b*b-c ;
        if (d >= 0)   // > 0
        {
           if(b>0){sd = c/( -b-std::sqrt(d));}
           else   {sd = -b + std::sqrt(d) ;}

          if ( sd >= 0 )   // > 0
          {
            if ( sd>dRmax ) // Avoid rounding errors due to precision issues on
            {               // 64 bits systems. Split long distance and recompute
              G4double fTerm = sd-std::fmod(sd,dRmax);
              sd = fTerm + DistanceToIn(p+fTerm*v,v);
            } 
            zi = p.z() + sd*v.z() ;

            if ( std::fabs(zi) <= tolODz )
            {
              if ( !fPhiFullCone )
              {
                xi     = p.x() + sd*v.x() ;
                yi     = p.y() + sd*v.y() ;
                ri     = rMinAv + zi*tanRMin ;
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri ;

                if (cosPsi >= cosHDPhiIT)
                { 
                  if ( sd > halfRadTolerance )  { snxt=sd; }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
                    Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin);
                    if ( Normal.dot(v) <= 0 )  { snxt = sd; }
                  } 
                }
              }
              else
              {
                if ( sd > halfRadTolerance )  { return sd; }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x() + sd*v.x() ;
                  yi     = p.y() + sd*v.y() ;
                  risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
                  Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin) ;
                  if ( Normal.dot(v) <= 0 )  { return sd; }
                }
              }
            }
          }
        }
      }
      else  if ( nt3 < -rin*kRadTolerance*secRMin )
      {
        // Within radius of inner cone (real or imaginary)
        // -> Try 2nd root, with checking intersection is with real cone
        // -> If check fails, try 1st root, also checking intersection is
        //    on real cone

        b = nt2/nt1 ;
        c = nt3/nt1 ;
        d = b*b - c ;

        if ( d >= 0 )  // > 0
        {
          if (b>0) { sd = c/(-b-std::sqrt(d)); }
          else     { sd = -b + std::sqrt(d);   }
          zi = p.z() + sd*v.z() ;
          ri = rMinAv + zi*tanRMin ;

          if ( ri > 0 )
          {
            if ( (sd >= 0) && (std::fabs(zi) <= tolODz) )  // sd > 0
            {
              if ( sd>dRmax ) // Avoid rounding errors due to precision issues
              {               // seen on 64 bits systems. Split and recompute
                G4double fTerm = sd-std::fmod(sd,dRmax);
                sd = fTerm + DistanceToIn(p+fTerm*v,v);
              } 
              if ( !fPhiFullCone )
              {
                xi     = p.x() + sd*v.x() ;
                yi     = p.y() + sd*v.y() ;
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri ;

                if (cosPsi >= cosHDPhiOT)
                {
                  if ( sd > halfRadTolerance )  { snxt=sd; }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
                    Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin);
                    if ( Normal.dot(v) <= 0 )  { snxt = sd; } 
                  }
                }
              }
              else
              {
                if( sd > halfRadTolerance )  { return sd; }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x() + sd*v.x() ;
                  yi     = p.y() + sd*v.y() ;
                  risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
                  Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin) ;
                  if ( Normal.dot(v) <= 0 )  { return sd; }
                } 
              }
            }
          }
          else
          {
            if (b>0) { sd = -b - std::sqrt(d);   }
            else     { sd = c/(-b+std::sqrt(d)); }
            zi = p.z() + sd*v.z() ;
            ri = rMinAv + zi*tanRMin ;

            if ( (sd >= 0) && (ri > 0) && (std::fabs(zi) <= tolODz) ) // sd>0
            {
              if ( sd>dRmax ) // Avoid rounding errors due to precision issues
              {               // seen on 64 bits systems. Split and recompute
                G4double fTerm = sd-std::fmod(sd,dRmax);
                sd = fTerm + DistanceToIn(p+fTerm*v,v);
              } 
              if ( !fPhiFullCone )
              {
                xi     = p.x() + sd*v.x() ;
                yi     = p.y() + sd*v.y() ;
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri ;

                if (cosPsi >= cosHDPhiIT)
                {
                  if ( sd > halfRadTolerance )  { snxt=sd; }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
                    Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin);
                    if ( Normal.dot(v) <= 0 )  { snxt = sd; } 
                  }
                }
              }
              else
              {
                if ( sd > halfRadTolerance )  { return sd; }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x() + sd*v.x() ;
                  yi     = p.y() + sd*v.y() ;
                  risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
                  Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin) ;
                  if ( Normal.dot(v) <= 0 )  { return sd; }
                } 
              }
            }
          }
        }
      }
      else
      {
        // Within kRadTol*0.5 of inner cone (real OR imaginary)
        // ----> Check not travelling through (=>0 to in)
        // ----> if not:
        //    -2nd root with validity check

        if ( std::fabs(p.z()) <= tolODz )
        {
          if ( nt2 > 0 )
          {
            // Inside inner real cone, heading outwards, inside z range

            if ( !fPhiFullCone )
            {
              cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/std::sqrt(t3) ;

              if (cosPsi >= cosHDPhiIT)  { return 0.0; }
            }
            else  { return 0.0; }
          }
          else
          {
            // Within z extent, but not travelling through
            // -> 2nd root or kInfinity if 1st root on imaginary cone

            b = nt2/nt1 ;
            c = nt3/nt1 ;
            d = b*b - c ;

            if ( d >= 0 )   // > 0
            {
              if (b>0) { sd = -b - std::sqrt(d);   }
              else     { sd = c/(-b+std::sqrt(d)); }
              zi = p.z() + sd*v.z() ;
              ri = rMinAv + zi*tanRMin ;
              
              if ( ri > 0 )   // 2nd root
              {
                if (b>0) { sd = c/(-b-std::sqrt(d)); }
                else     { sd = -b + std::sqrt(d);   }
                
                zi = p.z() + sd*v.z() ;

                if ( (sd >= 0) && (std::fabs(zi) <= tolODz) )  // sd>0
                {
                  if ( sd>dRmax ) // Avoid rounding errors due to precision issue
                  {               // seen on 64 bits systems. Split and recompute
                    G4double fTerm = sd-std::fmod(sd,dRmax);
                    sd = fTerm + DistanceToIn(p+fTerm*v,v);
                  } 
                  if ( !fPhiFullCone )
                  {
                    xi     = p.x() + sd*v.x() ;
                    yi     = p.y() + sd*v.y() ;
                    ri     = rMinAv + zi*tanRMin ;
                    cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri ;

                    if ( cosPsi >= cosHDPhiIT )  { snxt = sd; }
                  }
                  else  { return sd; }
                }
              }
              else  { return kInfinity; }
            }
          }
        }
        else   // 2nd root
        {
          b = nt2/nt1 ;
          c = nt3/nt1 ;
          d = b*b - c ;

          if ( d > 0 )
          {  
            if (b>0) { sd = c/(-b-std::sqrt(d)); }
            else     { sd = -b + std::sqrt(d) ;  }
            zi = p.z() + sd*v.z() ;

            if ( (sd >= 0) && (std::fabs(zi) <= tolODz) )  // sd>0
            {
              if ( sd>dRmax ) // Avoid rounding errors due to precision issues
              {               // seen on 64 bits systems. Split and recompute
                G4double fTerm = sd-std::fmod(sd,dRmax);
                sd = fTerm + DistanceToIn(p+fTerm*v,v);
              } 
              if ( !fPhiFullCone )
              {
                xi     = p.x() + sd*v.x();
                yi     = p.y() + sd*v.y();
                ri     = rMinAv + zi*tanRMin ;
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/ri;

                if (cosPsi >= cosHDPhiIT)  { snxt = sd; }
              }
              else  { return sd; }
            }
          }
        }
      }
    }
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to kCarTolerance*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> Should use some form of loop Construct

  if ( !fPhiFullCone )
  {
    // First phi surface (starting phi)

    Comp    = v.x()*sinSPhi - v.y()*cosSPhi ;
                    
    if ( Comp < 0 )    // Component in outwards normal dirn
    {
      Dist = (p.y()*cosSPhi - p.x()*sinSPhi) ;

      if (Dist < halfCarTolerance)
      {
        sd = Dist/Comp ;

        if ( sd < snxt )
        {
          if ( sd < 0 )  { sd = 0.0; }

          zi = p.z() + sd*v.z() ;

          if ( std::fabs(zi) <= tolODz )
          {
            xi        = p.x() + sd*v.x() ;
            yi        = p.y() + sd*v.y() ;
            rhoi2     = xi*xi + yi*yi ;
            tolORMin2 = (rMinOAv + zi*tanRMin)*(rMinOAv + zi*tanRMin) ;
            tolORMax2 = (rMaxOAv + zi*tanRMax)*(rMaxOAv + zi*tanRMax) ;

            if ( (rhoi2 >= tolORMin2) && (rhoi2 <= tolORMax2) )
            {
              // z and r intersections good - check intersecting with
              // correct half-plane

              if ((yi*cosCPhi - xi*sinCPhi) <= 0 )  { snxt = sd; }
            }
          }
        }
      }
    }

    // Second phi surface (Ending phi)

    Comp    = -(v.x()*sinEPhi - v.y()*cosEPhi) ;
        
    if ( Comp < 0 )   // Component in outwards normal dirn
    {
      Dist = -(p.y()*cosEPhi - p.x()*sinEPhi) ;
      if (Dist < halfCarTolerance)
      {
        sd = Dist/Comp ;

        if ( sd < snxt )
        {
          if ( sd < 0 )  { sd = 0.0; }

          zi = p.z() + sd*v.z() ;

          if (std::fabs(zi) <= tolODz)
          {
            xi        = p.x() + sd*v.x() ;
            yi        = p.y() + sd*v.y() ;
            rhoi2     = xi*xi + yi*yi ;
            tolORMin2 = (rMinOAv + zi*tanRMin)*(rMinOAv + zi*tanRMin) ;
            tolORMax2 = (rMaxOAv + zi*tanRMax)*(rMaxOAv + zi*tanRMax) ;

            if ( (rhoi2 >= tolORMin2) && (rhoi2 <= tolORMax2) )
            {
              // z and r intersections good - check intersecting with
              // correct half-plane

              if ( (yi*cosCPhi - xi*sinCPhi) >= 0.0 )  { snxt = sd; }
            }
          }
        }
      }
    }
  }
  if (snxt < halfCarTolerance)  { snxt = 0.; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////////////////
// 
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Cons::DistanceToIn(const G4ThreeVector& p) const
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ, safePhi, cosPsi ;
  G4double tanRMin, secRMin, pRMin ;
  G4double tanRMax, secRMax, pRMax ;

  rho   = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;
  safeZ = std::fabs(p.z()) - fDz ;

  if ( fRmin1 || fRmin2 )
  {
    tanRMin = (fRmin2 - fRmin1)*0.5/fDz ;
    secRMin = std::sqrt(1.0 + tanRMin*tanRMin) ;
    pRMin   = tanRMin*p.z() + (fRmin1 + fRmin2)*0.5 ;
    safeR1  = (pRMin - rho)/secRMin ;

    tanRMax = (fRmax2 - fRmax1)*0.5/fDz ;
    secRMax = std::sqrt(1.0 + tanRMax*tanRMax) ;
    pRMax   = tanRMax*p.z() + (fRmax1 + fRmax2)*0.5 ;
    safeR2  = (rho - pRMax)/secRMax ;

    if ( safeR1 > safeR2) { safe = safeR1; }
    else                  { safe = safeR2; }
  }
  else
  {
    tanRMax = (fRmax2 - fRmax1)*0.5/fDz ;
    secRMax = std::sqrt(1.0 + tanRMax*tanRMax) ;
    pRMax   = tanRMax*p.z() + (fRmax1 + fRmax2)*0.5 ;
    safe    = (rho - pRMax)/secRMax ;
  }
  if ( safeZ > safe )  { safe = safeZ; }

  if ( !fPhiFullCone && rho )
  {
    // Psi=angle from central phi to point

    cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/rho ;

    if ( cosPsi < std::cos(fDPhi*0.5) ) // Point lies outside phi range
    {
      if ( (p.y()*cosCPhi - p.x()*sinCPhi) <= 0.0 )
      {
        safePhi = std::fabs(p.x()*std::sin(fSPhi)-p.y()*std::cos(fSPhi));
      }
      else
      {
        safePhi = std::fabs(p.x()*sinEPhi-p.y()*cosEPhi);
      }
      if ( safePhi > safe )  { safe = safePhi; }
    }
  }
  if ( safe < 0.0 )  { safe = 0.0; }

  return safe ;
}

///////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from 'inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4Cons::DistanceToOut( const G4ThreeVector& p,
                                const G4ThreeVector& v,
                                const G4bool calcNorm,
                                      G4bool *validNorm,
                                      G4ThreeVector *n) const
{
  ESide side = kNull, sider = kNull, sidephi = kNull;

  G4double snxt,srd,sphi,pdist ;

  G4double tanRMax, secRMax, rMaxAv ;  // Data for outer cone
  G4double tanRMin, secRMin, rMinAv ;  // Data for inner cone

  G4double t1, t2, t3, rout, rin, nt1, nt2, nt3 ;
  G4double b, c, d, sr2, sr3 ;

  // Vars for intersection within tolerance

  ESide    sidetol = kNull ;
  G4double slentol = kInfinity ;

  // Vars for phi intersection:

  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, risec, vphi ;
  G4double zi, ri, deltaRoi2 ;

  // Z plane intersection

  if ( v.z() > 0.0 )
  {
    pdist = fDz - p.z() ;

    if (pdist > halfCarTolerance)
    {
      snxt = pdist/v.z() ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector(0,0,1) ;
        *validNorm = true ;
      }
      return  snxt = 0.0;
    }
  }
  else if ( v.z() < 0.0 )
  {
    pdist = fDz + p.z() ;

    if ( pdist > halfCarTolerance)
    {
      snxt = -pdist/v.z() ;
      side = kMZ ;
    }
    else
    {
      if ( calcNorm )
      {
        *n         = G4ThreeVector(0,0,-1) ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
  }
  else     // Travel perpendicular to z axis
  {
    snxt = kInfinity ;    
    side = kNull ;
  }

  // Radial Intersections
  //
  // Intersection with outer cone (possible return) and
  //                   inner cone (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=(a*z+b)^2
  //
  // where a=tanRMax or tanRMin
  //       b=rMaxAv  or rMinAv
  //
  // (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0 ;
  //     t1                        t2                      t3  
  //
  //  \--------u-------/       \-----------v----------/ \---------w--------/

  tanRMax = (fRmax2 - fRmax1)*0.5/fDz ;
  secRMax = std::sqrt(1.0 + tanRMax*tanRMax) ;
  rMaxAv  = (fRmax1 + fRmax2)*0.5 ;


  t1   = 1.0 - v.z()*v.z() ;      // since v normalised
  t2   = p.x()*v.x() + p.y()*v.y() ;
  t3   = p.x()*p.x() + p.y()*p.y() ;
  rout = tanRMax*p.z() + rMaxAv ;

  nt1 = t1 - (tanRMax*v.z())*(tanRMax*v.z()) ;
  nt2 = t2 - tanRMax*v.z()*rout ;
  nt3 = t3 - rout*rout ;

  if (v.z() > 0.0)
  {
    deltaRoi2 = snxt*snxt*t1 + 2*snxt*t2 + t3
                - fRmax2*(fRmax2 + kRadTolerance*secRMax);
  }
  else if ( v.z() < 0.0 )
  {
    deltaRoi2 = snxt*snxt*t1 + 2*snxt*t2 + t3
                - fRmax1*(fRmax1 + kRadTolerance*secRMax);
  }
  else
  {
    deltaRoi2 = 1.0;
  }

  if ( nt1 && (deltaRoi2 > 0.0) )  
  {
    // Equation quadratic => 2 roots : second root must be leaving

    b = nt2/nt1 ;
    c = nt3/nt1 ;
    d = b*b - c ;

    if ( d >= 0 )
    {
      // Check if on outer cone & heading outwards
      // NOTE: Should use rho-rout>-kRadTolerance*0.5
        
      if (nt3 > -halfRadTolerance && nt2 >= 0 )
      {
        if (calcNorm)
        {
          risec      = std::sqrt(t3)*secRMax ;
          *validNorm = true ;
          *n         = G4ThreeVector(p.x()/risec,p.y()/risec,-tanRMax/secRMax);
        }
        return snxt=0 ;
      }
      else
      {
        sider = kRMax  ;
        if (b>0) { srd = -b - std::sqrt(d);    }
        else     { srd = c/(-b+std::sqrt(d)) ; }

        zi    = p.z() + srd*v.z() ;
        ri    = tanRMax*zi + rMaxAv ;
          
        if ((ri >= 0) && (-halfRadTolerance <= srd) && (srd <= halfRadTolerance))
        {
          // An intersection within the tolerance
          //   we will Store it in case it is good -
          // 
          slentol = srd ;
          sidetol = kRMax ;
        }            
        if ( (ri < 0) || (srd < halfRadTolerance) )
        {
          // Safety: if both roots -ve ensure that srd cannot `win'
          //         distance to out

          if (b>0) { sr2 = c/(-b-std::sqrt(d)); }
          else     { sr2 = -b + std::sqrt(d);   }
          zi  = p.z() + sr2*v.z() ;
          ri  = tanRMax*zi + rMaxAv ;

          if ((ri >= 0) && (sr2 > halfRadTolerance))
          {
            srd = sr2;
          }
          else
          {
            srd = kInfinity ;

            if( (-halfRadTolerance <= sr2) && ( sr2 <= halfRadTolerance) )
            {
              // An intersection within the tolerance.
              // Storing it in case it is good.

              slentol = sr2 ;
              sidetol = kRMax ;
            }
          }
        }
      }
    }
    else
    {
      // No intersection with outer cone & not parallel
      // -> already outside, no intersection

      if ( calcNorm )
      {
        risec      = std::sqrt(t3)*secRMax;
        *validNorm = true;
        *n         = G4ThreeVector(p.x()/risec,p.y()/risec,-tanRMax/secRMax);
      }
      return snxt = 0.0 ;
    }
  }
  else if ( nt2 && (deltaRoi2 > 0.0) )
  {
    // Linear case (only one intersection) => point outside outer cone

    if ( calcNorm )
    {
      risec      = std::sqrt(t3)*secRMax;
      *validNorm = true;
      *n         = G4ThreeVector(p.x()/risec,p.y()/risec,-tanRMax/secRMax);
    }
    return snxt = 0.0 ;
  }
  else
  {
    // No intersection -> parallel to outer cone
    // => Z or inner cone intersection

    srd = kInfinity ;
  }

  // Check possible intersection within tolerance

  if ( slentol <= halfCarTolerance )
  {
    // An intersection within the tolerance was found.  
    // We must accept it only if the momentum points outwards.  
    //
    // G4ThreeVector ptTol ;  // The point of the intersection  
    // ptTol= p + slentol*v ;
    // ri=tanRMax*zi+rMaxAv ;
    //
    // Calculate a normal vector,  as below

    xi    = p.x() + slentol*v.x();
    yi    = p.y() + slentol*v.y();
    risec = std::sqrt(xi*xi + yi*yi)*secRMax;
    G4ThreeVector Normal = G4ThreeVector(xi/risec,yi/risec,-tanRMax/secRMax);

    if ( Normal.dot(v) > 0 )    // We will leave the Cone immediatelly
    {
      if ( calcNorm ) 
      {
        *n         = Normal.unit() ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
    else // On the surface, but not heading out so we ignore this intersection
    {    //                                        (as it is within tolerance).
      slentol = kInfinity ;
    }
  }

  // Inner Cone intersection

  if ( fRmin1 || fRmin2 )
  {
    tanRMin = (fRmin2 - fRmin1)*0.5/fDz ;
    nt1     = t1 - (tanRMin*v.z())*(tanRMin*v.z()) ;

    if ( nt1 )
    {
      secRMin = std::sqrt(1.0 + tanRMin*tanRMin) ;
      rMinAv  = (fRmin1 + fRmin2)*0.5 ;    
      rin     = tanRMin*p.z() + rMinAv ;
      nt2     = t2 - tanRMin*v.z()*rin ;
      nt3     = t3 - rin*rin ;
      
      // Equation quadratic => 2 roots : first root must be leaving

      b = nt2/nt1 ;
      c = nt3/nt1 ;
      d = b*b - c ;

      if ( d >= 0.0 )
      {
        // NOTE: should be rho-rin<kRadTolerance*0.5,
        //       but using squared versions for efficiency

        if (nt3 < kRadTolerance*(rin + kRadTolerance*0.25)) 
        {
          if ( nt2 < 0.0 )
          {
            if (calcNorm)  { *validNorm = false; }
            return          snxt      = 0.0;
          }
        }
        else
        {
          if (b>0) { sr2 = -b - std::sqrt(d);   }
          else     { sr2 = c/(-b+std::sqrt(d)); }
          zi  = p.z() + sr2*v.z() ;
          ri  = tanRMin*zi + rMinAv ;

          if( (ri>=0.0)&&(-halfRadTolerance<=sr2)&&(sr2<=halfRadTolerance) )
          {
            // An intersection within the tolerance
            // storing it in case it is good.

            slentol = sr2 ;
            sidetol = kRMax ;
          }
          if( (ri<0) || (sr2 < halfRadTolerance) )
          {
            if (b>0) { sr3 = c/(-b-std::sqrt(d)); }
            else     { sr3 = -b + std::sqrt(d) ;  }

            // Safety: if both roots -ve ensure that srd cannot `win'
            //         distancetoout

            if  ( sr3 > halfRadTolerance )
            {
              if( sr3 < srd )
              {
                zi = p.z() + sr3*v.z() ;
                ri = tanRMin*zi + rMinAv ;

                if ( ri >= 0.0 )
                {
                  srd=sr3 ;
                  sider=kRMin ;
                }
              } 
            }
            else if ( sr3 > -halfRadTolerance )
            {
              // Intersection in tolerance. Store to check if it's good

              slentol = sr3 ;
              sidetol = kRMin ;
            }
          }
          else if ( (sr2 < srd) && (sr2 > halfCarTolerance) )
          {
            srd   = sr2 ;
            sider = kRMin ;
          }
          else if (sr2 > -halfCarTolerance)
          {
            // Intersection in tolerance. Store to check if it's good

            slentol = sr2 ;
            sidetol = kRMin ;
          }    
          if( slentol <= halfCarTolerance  )
          {
            // An intersection within the tolerance was found. 
            // We must accept it only if  the momentum points outwards. 

            G4ThreeVector Normal ; 
            
            // Calculate a normal vector,  as below

            xi     = p.x() + slentol*v.x() ;
            yi     = p.y() + slentol*v.y() ;
            if( sidetol==kRMax )
            {
              risec  = std::sqrt(xi*xi + yi*yi)*secRMax ;
              Normal = G4ThreeVector(xi/risec,yi/risec,-tanRMax/secRMax) ;
            }
            else
            {
              risec  = std::sqrt(xi*xi + yi*yi)*secRMin ;
              Normal = G4ThreeVector(-xi/risec,-yi/risec,tanRMin/secRMin) ;
            }
            if( Normal.dot(v) > 0 )
            {
              // We will leave the cone immediately

              if( calcNorm ) 
              {
                *n         = Normal.unit() ;
                *validNorm = true ;
              }
              return snxt = 0.0 ;
            }
            else 
            { 
              // On the surface, but not heading out so we ignore this
              // intersection (as it is within tolerance). 

              slentol = kInfinity ;
            }        
          }
        }
      }
    }
  }

  // Linear case => point outside inner cone ---> outer cone intersect
  //
  // Phi Intersection
  
  if ( !fPhiFullCone )
  {
    // add angle calculation with correction 
    // of the difference in domain of atan2 and Sphi

    vphi = std::atan2(v.y(),v.x()) ;

    if ( vphi < fSPhi - halfAngTolerance  )              { vphi += twopi; }
    else if ( vphi > fSPhi + fDPhi + halfAngTolerance )  { vphi -= twopi; }

    if ( p.x() || p.y() )   // Check if on z axis (rho not needed later)
    {
      // pDist -ve when inside

      pDistS = p.x()*sinSPhi - p.y()*cosSPhi ;
      pDistE = -p.x()*sinEPhi + p.y()*cosEPhi ;

      // Comp -ve when in direction of outwards normal

      compS = -sinSPhi*v.x() + cosSPhi*v.y() ;
      compE = sinEPhi*v.x() - cosEPhi*v.y() ;

      sidephi = kNull ;

      if( ( (fDPhi <= pi) && ( (pDistS <= halfCarTolerance)
                            && (pDistE <= halfCarTolerance) ) )
         || ( (fDPhi >  pi) && !((pDistS >  halfCarTolerance)
                              && (pDistE >  halfCarTolerance) ) )  )
      {
        // Inside both phi *full* planes
        if ( compS < 0 )
        {
          sphi = pDistS/compS ;
          if (sphi >= -halfCarTolerance)
          {
            xi = p.x() + sphi*v.x() ;
            yi = p.y() + sphi*v.y() ;

            // Check intersecting with correct half-plane
            // (if not -> no intersect)
            //
            if ( (std::fabs(xi)<=kCarTolerance)
              && (std::fabs(yi)<=kCarTolerance) )
            {
              sidephi= kSPhi;
              if ( ( fSPhi-halfAngTolerance <= vphi )
                && ( fSPhi+fDPhi+halfAngTolerance >=vphi ) )
              {
                sphi = kInfinity;
              }
            }
            else
            if ( (yi*cosCPhi-xi*sinCPhi)>=0 )
            {
              sphi = kInfinity ;
            }
            else
            {
              sidephi = kSPhi ;
              if ( pDistS > -halfCarTolerance )
              {
                sphi = 0.0 ; // Leave by sphi immediately
              }    
            }       
          }
          else
          {
            sphi = kInfinity ;
          }
        }
        else
        {
          sphi = kInfinity ;
        }

        if ( compE < 0 )
        {
          sphi2 = pDistE/compE ;

          // Only check further if < starting phi intersection
          //
          if ( (sphi2 > -halfCarTolerance) && (sphi2 < sphi) )
          {
            xi = p.x() + sphi2*v.x() ;
            yi = p.y() + sphi2*v.y() ;

            // Check intersecting with correct half-plane

            if ( (std::fabs(xi)<=kCarTolerance)
              && (std::fabs(yi)<=kCarTolerance) )
            {
              // Leaving via ending phi

              if(!( (fSPhi-halfAngTolerance <= vphi)
                 && (fSPhi+fDPhi+halfAngTolerance >= vphi) ) )
              {
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
                else                                { sphi = 0.0; }
              }
            }
            else // Check intersecting with correct half-plane
            if ( yi*cosCPhi-xi*sinCPhi >= 0 )
            {
              // Leaving via ending phi

              sidephi = kEPhi ;
              if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
              else                                { sphi = 0.0; }
            }
          }
        }
      }
      else
      {
        sphi = kInfinity ;
      }
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if ( (fSPhi-halfAngTolerance <= vphi)
        && (vphi <= fSPhi+fDPhi+halfAngTolerance) )
      {
        sphi = kInfinity ;
      }
      else
      {
        sidephi = kSPhi  ;   // arbitrary 
        sphi    = 0.0 ;
      }
    }      
    if ( sphi < snxt )  // Order intersecttions
    {
      snxt=sphi ;
      side=sidephi ;
    }
  }
  if ( srd < snxt )  // Order intersections
  {
    snxt = srd   ;
    side = sider ;
  }
  if (calcNorm)
  {
    switch(side)
    {                     // Note: returned vector not normalised
      case kRMax:         // (divide by frmax for unit vector)
        xi         = p.x() + snxt*v.x() ;
        yi         = p.y() + snxt*v.y() ;
        risec      = std::sqrt(xi*xi + yi*yi)*secRMax ;
        *n         = G4ThreeVector(xi/risec,yi/risec,-tanRMax/secRMax) ;
        *validNorm = true ;
        break ;
      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;
      case kSPhi:
        if ( fDPhi <= pi )
        {
          *n         = G4ThreeVector(sinSPhi, -cosSPhi, 0);
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;
      case kEPhi:
        if ( fDPhi <= pi )
        {
          *n = G4ThreeVector(-sinEPhi, cosEPhi, 0);
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;
      case kPZ:
        *n         = G4ThreeVector(0,0,1) ;
        *validNorm = true ;
        break ;
      case kMZ:
        *n         = G4ThreeVector(0,0,-1) ;
        *validNorm = true ;
        break ;
      default:
        G4cout << G4endl ;
        DumpInfo();
        std::ostringstream message;
        G4int oldprc = message.precision(16) ;
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x() = "   << p.x()/mm << " mm" << G4endl
                << "p.y() = "   << p.y()/mm << " mm" << G4endl
                << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
                << "pho at z = "   << std::sqrt( p.x()*p.x()+p.y()*p.y() )/mm
                << " mm" << G4endl << G4endl ;
        if( p.x() != 0. || p.y() != 0.)
        {
           message << "point phi = "   << std::atan2(p.y(),p.x())/degree
                   << " degree" << G4endl << G4endl ; 
        }
        message << "Direction:" << G4endl << G4endl
                << "v.x() = "   << v.x() << G4endl
                << "v.y() = "   << v.y() << G4endl
                << "v.z() = "   << v.z() << G4endl<< G4endl
                << "Proposed distance :" << G4endl<< G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl ;
        message.precision(oldprc) ;
        G4Exception("G4Cons::DistanceToOut(p,v,..)","GeomSolids1002",
                    JustWarning, message) ;
        break ;
    }
  }
  if (snxt < halfCarTolerance)  { snxt = 0.; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Cons::DistanceToOut(const G4ThreeVector& p) const
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ, safePhi;
  G4double tanRMin, secRMin, pRMin;
  G4double tanRMax, secRMax, pRMax;

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
    G4int oldprc=G4cout.precision(16) ;
    G4cout << G4endl ;
    DumpInfo();
    G4cout << "Position:"  << G4endl << G4endl ;
    G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
    G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
    G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
    G4cout << "pho at z = "   << std::sqrt( p.x()*p.x()+p.y()*p.y() )/mm
           << " mm" << G4endl << G4endl ;
    if( (p.x() != 0.) || (p.x() != 0.) )
    {
      G4cout << "point phi = "   << std::atan2(p.y(),p.x())/degree
             << " degree" << G4endl << G4endl ; 
    }
    G4cout.precision(oldprc) ;
    G4Exception("G4Cons::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, "Point p is outside !?" );
  }
#endif

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;
  safeZ = fDz - std::fabs(p.z()) ;

  if (fRmin1 || fRmin2)
  {
    tanRMin = (fRmin2 - fRmin1)*0.5/fDz ;
    secRMin = std::sqrt(1.0 + tanRMin*tanRMin) ;
    pRMin   = tanRMin*p.z() + (fRmin1 + fRmin2)*0.5 ;
    safeR1  = (rho - pRMin)/secRMin ;
  }
  else
  {
    safeR1 = kInfinity ;
  }

  tanRMax = (fRmax2 - fRmax1)*0.5/fDz ;
  secRMax = std::sqrt(1.0 + tanRMax*tanRMax) ;
  pRMax   = tanRMax*p.z() + (fRmax1+fRmax2)*0.5 ;
  safeR2  = (pRMax - rho)/secRMax ;

  if (safeR1 < safeR2)  { safe = safeR1; }
  else                  { safe = safeR2; }
  if (safeZ < safe)     { safe = safeZ ; }

  // Check if phi divided, Calc distances closest phi plane

  if (!fPhiFullCone)
  {
    // Above/below central phi of G4Cons?

    if ( (p.y()*cosCPhi - p.x()*sinCPhi) <= 0 )
    {
      safePhi = -(p.x()*sinSPhi - p.y()*cosSPhi) ;
    }
    else
    {
      safePhi = (p.x()*sinEPhi - p.y()*cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi; }
  }
  if ( safe < 0 )  { safe = 0; }

  return safe ;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Cons::GetEntityType() const
{
  return G4String("G4Cons");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Cons::Clone() const
{
  return new G4Cons(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Cons::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Cons\n"
     << " Parameters: \n"
     << "   inside  -fDz radius: "  << fRmin1/mm << " mm \n"
     << "   outside -fDz radius: "  << fRmax1/mm << " mm \n"
     << "   inside  +fDz radius: "  << fRmin2/mm << " mm \n"
     << "   outside +fDz radius: "  << fRmax2/mm << " mm \n"
     << "   half length in Z   : "  << fDz/mm << " mm \n"
     << "   starting angle of segment: " << fSPhi/degree << " degrees \n"
     << "   delta angle of segment   : " << fDPhi/degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}



/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Cons::GetPointOnSurface() const
{   
  // declare working variables
  //
  G4double Aone, Atwo, Athree, Afour, Afive, slin, slout, phi;
  G4double zRand, cosu, sinu, rRand1, rRand2, chose, rone, rtwo, qone, qtwo;
  rone = (fRmax1-fRmax2)/(2.*fDz);
  rtwo = (fRmin1-fRmin2)/(2.*fDz);
  qone=0.; qtwo=0.;
  if(fRmax1!=fRmax2) { qone = fDz*(fRmax1+fRmax2)/(fRmax1-fRmax2); }
  if(fRmin1!=fRmin2) { qtwo = fDz*(fRmin1+fRmin2)/(fRmin1-fRmin2); }
  slin   = std::sqrt(sqr(fRmin1-fRmin2)+sqr(2.*fDz));
  slout  = std::sqrt(sqr(fRmax1-fRmax2)+sqr(2.*fDz));
  Aone   = 0.5*fDPhi*(fRmax2 + fRmax1)*slout;       
  Atwo   = 0.5*fDPhi*(fRmin2 + fRmin1)*slin;
  Athree = 0.5*fDPhi*(fRmax1*fRmax1-fRmin1*fRmin1); 
  Afour  = 0.5*fDPhi*(fRmax2*fRmax2-fRmin2*fRmin2);
  Afive  = fDz*(fRmax1-fRmin1+fRmax2-fRmin2);
  
  phi    = G4RandFlat::shoot(fSPhi,fSPhi+fDPhi);
  cosu   = std::cos(phi);  sinu = std::sin(phi);
  rRand1 = GetRadiusInRing(fRmin1, fRmax1);
  rRand2 = GetRadiusInRing(fRmin2, fRmax2);
  
  if ( (fSPhi == 0.) && fPhiFullCone )  { Afive = 0.; }
  chose  = G4RandFlat::shoot(0.,Aone+Atwo+Athree+Afour+2.*Afive);
 
  if( (chose >= 0.) && (chose < Aone) )
  {
    if(fRmin1 != fRmin2)
    {
      zRand = G4RandFlat::shoot(-1.*fDz,fDz); 
      return G4ThreeVector (rtwo*cosu*(qtwo-zRand),
                            rtwo*sinu*(qtwo-zRand), zRand);
    }
    else
    {
      return G4ThreeVector(fRmin1*cosu, fRmin2*sinu,
                           G4RandFlat::shoot(-1.*fDz,fDz));
    }
  }
  else if( (chose >= Aone) && (chose <= Aone + Atwo) )
  {
    if(fRmax1 != fRmax2)
    {
      zRand = G4RandFlat::shoot(-1.*fDz,fDz); 
      return G4ThreeVector (rone*cosu*(qone-zRand),
                            rone*sinu*(qone-zRand), zRand);
    }    
    else
    {
      return G4ThreeVector(fRmax1*cosu, fRmax2*sinu,
                           G4RandFlat::shoot(-1.*fDz,fDz));
    }
  }
  else if( (chose >= Aone + Atwo) && (chose < Aone + Atwo + Athree) )
  {
    return G4ThreeVector (rRand1*cosu, rRand1*sinu, -1*fDz);
  }
  else if( (chose >= Aone + Atwo + Athree)
        && (chose < Aone + Atwo + Athree + Afour) )
  {
    return G4ThreeVector (rRand2*cosu,rRand2*sinu,fDz);
  }
  else if( (chose >= Aone + Atwo + Athree + Afour)
        && (chose < Aone + Atwo + Athree + Afour + Afive) )
  {
    zRand  = G4RandFlat::shoot(-1.*fDz,fDz);
    rRand1 = G4RandFlat::shoot(fRmin2-((zRand-fDz)/(2.*fDz))*(fRmin1-fRmin2),
                               fRmax2-((zRand-fDz)/(2.*fDz))*(fRmax1-fRmax2)); 
    return G4ThreeVector (rRand1*std::cos(fSPhi),
                          rRand1*std::sin(fSPhi), zRand);
  }
  else
  { 
    zRand  = G4RandFlat::shoot(-1.*fDz,fDz);
    rRand1 = G4RandFlat::shoot(fRmin2-((zRand-fDz)/(2.*fDz))*(fRmin1-fRmin2),
                               fRmax2-((zRand-fDz)/(2.*fDz))*(fRmax1-fRmax2)); 
    return G4ThreeVector (rRand1*std::cos(fSPhi+fDPhi),
                          rRand1*std::sin(fSPhi+fDPhi), zRand);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Cons::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Cons::CreatePolyhedron () const
{
  return new G4PolyhedronCons(fRmin1,fRmax1,fRmin2,fRmax2,fDz,fSPhi,fDPhi);
}

#endif
