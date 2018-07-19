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
// $Id: G4Sphere.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// class G4Sphere
//
// Implementation for G4Sphere class
//
// History:
//
// 26.10.16 E.Tcherniaev: re-implemented CalculateExtent() using
//                      G4BoundingEnvelope, removed CreateRotatedVertices()
// 05.04.12 M.Kelsey:   GetPointOnSurface() throw flat in cos(theta), sqrt(r)
// 14.09.09 T.Nikitina: fix for phi section in DistanceToOut(p,v,..),as for
//                      G4Tubs,G4Cons 
// 26.03.09 G.Cosmo   : optimisations and uniform use of local radial tolerance
// 12.06.08 V.Grichine: fix for theta intersections in DistanceToOut(p,v,...)
// 22.07.05 O.Link    : Added check for intersection with double cone
// 03.05.05 V.Grichine: SurfaceNormal(p) according to J. Apostolakis proposal
// 16.09.04 V.Grichine: bug fixed in SurfaceNormal(p), theta normals
// 16.07.04 V.Grichine: bug fixed in DistanceToOut(p,v), Rmin go outside
// 02.06.04 V.Grichine: bug fixed in DistanceToIn(p,v), on Rmax,Rmin go inside
// 30.10.03 J.Apostolakis: new algorithm in Inside for SPhi-sections
// 29.10.03 J.Apostolakis: fix in Inside for SPhi-0.5*kAngTol < phi<SPhi, SPhi<0
// 19.06.02 V.Grichine: bug fixed in Inside(p), && -> && fDTheta - kAngTolerance
// 30.01.02 V.Grichine: bug fixed in Inside(p), && -> || at l.451
// 06.03.00 V.Grichine: modifications in Distance ToOut(p,v,...)
// 18.11.99 V.Grichine: side = kNull in Distance ToOut(p,v,...)
// 25.11.98 V.Grichine: bug fixed in DistanceToIn(p,v), phi intersections
// 12.11.98 V.Grichine: bug fixed in DistanceToIn(p,v), theta intersections
// 09.10.98 V.Grichine: modifications in DistanceToOut(p,v,...)
// 17.09.96 V.Grichine: final modifications to commit
// 28.03.94 P.Kent: old C++ code converted to tolerant geometry
// --------------------------------------------------------------------

#include "G4Sphere.hh"

#if !defined(G4GEOM_USE_USPHERE)

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

using namespace CLHEP;

// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kSTheta,kETheta};

// used by normal

enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNSTheta,kNETheta};

////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4Sphere::G4Sphere( const G4String& pName,
                          G4double pRmin, G4double pRmax,
                          G4double pSPhi, G4double pDPhi,
                          G4double pSTheta, G4double pDTheta )
  : G4CSGSolid(pName), fEpsilon(2.e-11), fSPhi(0.0),
    fFullPhiSphere(true), fFullThetaSphere(true)
{
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();

  halfCarTolerance = 0.5*kCarTolerance;
  halfAngTolerance = 0.5*kAngTolerance;

  // Check radii and set radial tolerances

  if ( (pRmin >= pRmax) || (pRmax < 1.1*kRadTolerance) || (pRmin < 0) )
  {
    std::ostringstream message;
    message << "Invalid radii for Solid: " << GetName() << G4endl
            << "        pRmin = " << pRmin << ", pRmax = " << pRmax;
    G4Exception("G4Sphere::G4Sphere()", "GeomSolids0002",
                FatalException, message);
  }
  fRmin=pRmin; fRmax=pRmax;
  fRminTolerance = (fRmin) ? std::max( kRadTolerance, fEpsilon*fRmin ) : 0;
  fRmaxTolerance = std::max( kRadTolerance, fEpsilon*fRmax );

  // Check angles

  CheckPhiAngles(pSPhi, pDPhi);
  CheckThetaAngles(pSTheta, pDTheta);
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Sphere::G4Sphere( __void__& a )
  : G4CSGSolid(a), fRminTolerance(0.), fRmaxTolerance(0.),
    kAngTolerance(0.), kRadTolerance(0.), fEpsilon(0.),
    fRmin(0.), fRmax(0.), fSPhi(0.), fDPhi(0.), fSTheta(0.),
    fDTheta(0.), sinCPhi(0.), cosCPhi(0.), cosHDPhiOT(0.), cosHDPhiIT(0.),
    sinSPhi(0.), cosSPhi(0.), sinEPhi(0.), cosEPhi(0.), hDPhi(0.), cPhi(0.),
    ePhi(0.), sinSTheta(0.), cosSTheta(0.), sinETheta(0.), cosETheta(0.),
    tanSTheta(0.), tanSTheta2(0.), tanETheta(0.), tanETheta2(0.), eTheta(0.),
    fFullPhiSphere(false), fFullThetaSphere(false), fFullSphere(true),
    halfCarTolerance(0.), halfAngTolerance(0.)
{
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4Sphere::~G4Sphere()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Sphere::G4Sphere(const G4Sphere& rhs)
  : G4CSGSolid(rhs), fRminTolerance(rhs.fRminTolerance),
    fRmaxTolerance(rhs.fRmaxTolerance), kAngTolerance(rhs.kAngTolerance),
    kRadTolerance(rhs.kRadTolerance), fEpsilon(rhs.fEpsilon),
    fRmin(rhs.fRmin), fRmax(rhs.fRmax), fSPhi(rhs.fSPhi), fDPhi(rhs.fDPhi),
    fSTheta(rhs.fSTheta), fDTheta(rhs.fDTheta),
    sinCPhi(rhs.sinCPhi), cosCPhi(rhs.cosCPhi),
    cosHDPhiOT(rhs.cosHDPhiOT), cosHDPhiIT(rhs.cosHDPhiIT),
    sinSPhi(rhs.sinSPhi), cosSPhi(rhs.cosSPhi),
    sinEPhi(rhs.sinEPhi), cosEPhi(rhs.cosEPhi),
    hDPhi(rhs.hDPhi), cPhi(rhs.cPhi), ePhi(rhs.ePhi),
    sinSTheta(rhs.sinSTheta), cosSTheta(rhs.cosSTheta),
    sinETheta(rhs.sinETheta), cosETheta(rhs.cosETheta),
    tanSTheta(rhs.tanSTheta), tanSTheta2(rhs.tanSTheta2),
    tanETheta(rhs.tanETheta), tanETheta2(rhs.tanETheta2), eTheta(rhs.eTheta),
    fFullPhiSphere(rhs.fFullPhiSphere), fFullThetaSphere(rhs.fFullThetaSphere),
    fFullSphere(rhs.fFullSphere),
    halfCarTolerance(rhs.halfCarTolerance),
    halfAngTolerance(rhs.halfAngTolerance)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Sphere& G4Sphere::operator = (const G4Sphere& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fRminTolerance = rhs.fRminTolerance; fRmaxTolerance = rhs.fRmaxTolerance;
   kAngTolerance = rhs.kAngTolerance; kRadTolerance = rhs.kRadTolerance;
   fEpsilon = rhs.fEpsilon; fRmin = rhs.fRmin; fRmax = rhs.fRmax;
   fSPhi = rhs.fSPhi; fDPhi = rhs.fDPhi; fSTheta = rhs.fSTheta;
   fDTheta = rhs.fDTheta; sinCPhi = rhs.sinCPhi; cosCPhi = rhs.cosCPhi;
   cosHDPhiOT = rhs.cosHDPhiOT; cosHDPhiIT = rhs.cosHDPhiIT;
   sinSPhi = rhs.sinSPhi; cosSPhi = rhs.cosSPhi;
   sinEPhi = rhs.sinEPhi; cosEPhi = rhs.cosEPhi;
   hDPhi = rhs.hDPhi; cPhi = rhs.cPhi; ePhi = rhs.ePhi;
   sinSTheta = rhs.sinSTheta; cosSTheta = rhs.cosSTheta;
   sinETheta = rhs.sinETheta; cosETheta = rhs.cosETheta;
   tanSTheta = rhs.tanSTheta; tanSTheta2 = rhs.tanSTheta2;
   tanETheta = rhs.tanETheta; tanETheta2 = rhs.tanETheta2;
   eTheta = rhs.eTheta; fFullPhiSphere = rhs.fFullPhiSphere;
   fFullThetaSphere = rhs.fFullThetaSphere; fFullSphere = rhs.fFullSphere;
   halfCarTolerance = rhs.halfCarTolerance;
   halfAngTolerance = rhs.halfAngTolerance;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Sphere::ComputeDimensions(       G4VPVParameterisation* p,
                                  const G4int n,
                                  const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Sphere::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
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
    G4Exception("G4Sphere::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Sphere::CalculateExtent( const EAxis pAxis,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax ) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
// Split into radius, phi, theta checks
// Each check modifies 'in', or returns as approprate

EInside G4Sphere::Inside( const G4ThreeVector& p ) const
{
  G4double rho,rho2,rad2,tolRMin,tolRMax;
  G4double pPhi,pTheta;
  EInside in = kOutside;

  const G4double halfRmaxTolerance = fRmaxTolerance*0.5;
  const G4double halfRminTolerance = fRminTolerance*0.5;
  const G4double Rmax_minus = fRmax - halfRmaxTolerance;
  const G4double Rmin_plus  = (fRmin > 0) ? fRmin+halfRminTolerance : 0;

  rho2 = p.x()*p.x() + p.y()*p.y() ;
  rad2 = rho2 + p.z()*p.z() ;

  // Check radial surfaces. Sets 'in'

  tolRMin = Rmin_plus;
  tolRMax = Rmax_minus;

  if(rad2 == 0.0)
  { 
    if (fRmin > 0.0)
    {
      return in = kOutside;
    }
    if ( (!fFullPhiSphere) || (!fFullThetaSphere) )
    {
      return in = kSurface;
    }
    else
    {
      return in = kInside; 
    }
  }

  if ( (rad2 <= Rmax_minus*Rmax_minus) && (rad2 >= Rmin_plus*Rmin_plus) )
  {
    in = kInside;
  }
  else
  {
    tolRMax = fRmax + halfRmaxTolerance;                  // outside case
    tolRMin = std::max(fRmin-halfRminTolerance, 0.);      // outside case
    if ( (rad2 <= tolRMax*tolRMax) && (rad2 >= tolRMin*tolRMin) )
    {
      in = kSurface;
    }
    else
    {
      return in = kOutside;
    }
  }

  // Phi boundaries   : Do not check if it has no phi boundary!

  if ( !fFullPhiSphere && rho2 )  // [fDPhi < twopi] and [p.x or p.y]
  {
    pPhi = std::atan2(p.y(),p.x()) ;

    if      ( pPhi < fSPhi - halfAngTolerance  ) { pPhi += twopi; }
    else if ( pPhi > ePhi + halfAngTolerance )   { pPhi -= twopi; }
    
    if ( (pPhi < fSPhi - halfAngTolerance)
      || (pPhi > ePhi + halfAngTolerance) )      { return in = kOutside; }
    
    else if (in == kInside)  // else it's kSurface anyway already
    {
      if ( (pPhi < fSPhi + halfAngTolerance)
        || (pPhi > ePhi - halfAngTolerance) )    { in = kSurface; }     
    }
  }

  // Theta bondaries
  
  if ( (rho2 || p.z()) && (!fFullThetaSphere) )
  {
    rho    = std::sqrt(rho2);
    pTheta = std::atan2(rho,p.z());

    if ( in == kInside )
    {
      if ( ((fSTheta > 0.0) && (pTheta < fSTheta + halfAngTolerance))
        || ((eTheta < pi) && (pTheta > eTheta - halfAngTolerance)) )
      {
        if ( (( (fSTheta>0.0)&&(pTheta>=fSTheta-halfAngTolerance) )
             || (fSTheta == 0.0) )
          && ((eTheta==pi)||(pTheta <= eTheta + halfAngTolerance) ) )
        {
          in = kSurface;
        }
        else
        {
          in = kOutside;
        }
      }
    }
    else
    {
        if ( ((fSTheta > 0.0)&&(pTheta < fSTheta - halfAngTolerance))
           ||((eTheta < pi  )&&(pTheta > eTheta + halfAngTolerance)) )
      {
        in = kOutside;
      }
    }
  }
  return in;
}

/////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Sphere::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4int noSurfaces = 0;  
  G4double rho, rho2, radius, pTheta, pPhi=0.;
  G4double distRMin = kInfinity;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;
  G4double distSTheta = kInfinity, distETheta = kInfinity;
  G4ThreeVector nR, nPs, nPe, nTs, nTe, nZ(0.,0.,1.);
  G4ThreeVector norm, sumnorm(0.,0.,0.);

  rho2 = p.x()*p.x()+p.y()*p.y();
  radius = std::sqrt(rho2+p.z()*p.z());
  rho  = std::sqrt(rho2);

  G4double    distRMax = std::fabs(radius-fRmax);
  if (fRmin)  distRMin = std::fabs(radius-fRmin);
    
  if ( rho && !fFullSphere )
  {
    pPhi = std::atan2(p.y(),p.x());

    if (pPhi < fSPhi-halfAngTolerance)     { pPhi += twopi; }
    else if (pPhi > ePhi+halfAngTolerance) { pPhi -= twopi; }
  }
  if ( !fFullPhiSphere )
  {
    if ( rho )
    {
      distSPhi = std::fabs( pPhi-fSPhi ); 
      distEPhi = std::fabs( pPhi-ePhi ); 
    }
    else if( !fRmin )
    {
      distSPhi = 0.; 
      distEPhi = 0.; 
    }
    nPs = G4ThreeVector(sinSPhi,-cosSPhi,0);
    nPe = G4ThreeVector(-sinEPhi,cosEPhi,0);
  }        
  if ( !fFullThetaSphere )
  {
    if ( rho )
    {
      pTheta     = std::atan2(rho,p.z());
      distSTheta = std::fabs(pTheta-fSTheta); 
      distETheta = std::fabs(pTheta-eTheta);
 
      nTs = G4ThreeVector(-cosSTheta*p.x()/rho,
                          -cosSTheta*p.y()/rho,
                           sinSTheta          );

      nTe = G4ThreeVector( cosETheta*p.x()/rho,
                           cosETheta*p.y()/rho,
                          -sinETheta          );    
    }
    else if( !fRmin )
    {
      if ( fSTheta )  
      {              
        distSTheta = 0.;
        nTs = G4ThreeVector(0.,0.,-1.);
      }
      if ( eTheta < pi )
      {              
        distETheta = 0.;
        nTe = G4ThreeVector(0.,0.,1.);
      }
    }    
  }
  if( radius )  { nR = G4ThreeVector(p.x()/radius,p.y()/radius,p.z()/radius); }

  if( distRMax <= halfCarTolerance )
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  if( fRmin && (distRMin <= halfCarTolerance) )
  {
    noSurfaces ++;
    sumnorm -= nR;
  }
  if( !fFullPhiSphere )   
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
  if ( !fFullThetaSphere )
  {
    if ((distSTheta <= halfAngTolerance) && (fSTheta > 0.))
    {
      noSurfaces ++;
      if ((radius <= halfCarTolerance) && fFullPhiSphere)  { sumnorm += nZ;  }
      else                                                 { sumnorm += nTs; }
    }
    if ((distETheta <= halfAngTolerance) && (eTheta < pi)) 
    {
      noSurfaces ++;
      if ((radius <= halfCarTolerance) && fFullPhiSphere)  { sumnorm -= nZ;  }
      else                                                 { sumnorm += nTe; }
      if(sumnorm.z() == 0.)  { sumnorm += nZ; }
    }
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4Sphere::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" ); 
#endif
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = sumnorm.unit(); }
  return norm;
}


/////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Sphere::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  ENorm side;
  G4ThreeVector norm;
  G4double rho,rho2,radius,pPhi,pTheta;
  G4double distRMin,distRMax,distSPhi,distEPhi,
           distSTheta,distETheta,distMin;

  rho2=p.x()*p.x()+p.y()*p.y();
  radius=std::sqrt(rho2+p.z()*p.z());
  rho=std::sqrt(rho2);

  //
  // Distance to r shells
  //

  distRMax=std::fabs(radius-fRmax);
  if (fRmin)
  {
    distRMin=std::fabs(radius-fRmin);
      
    if (distRMin<distRMax)
    {
      distMin=distRMin;
      side=kNRMin;
    }
    else
    {   
      distMin=distRMax;
      side=kNRMax;
    }
  }
  else
  {
    distMin=distRMax;
    side=kNRMax;
  }

  //
  // Distance to phi planes
  //
  // Protected against (0,0,z) 
    
  pPhi = std::atan2(p.y(),p.x());
  if (pPhi<0) { pPhi += twopi; }

  if (!fFullPhiSphere && rho)
  {
    if (fSPhi<0)
    {
      distSPhi=std::fabs(pPhi-(fSPhi+twopi))*rho;
    }
    else
    {
      distSPhi=std::fabs(pPhi-fSPhi)*rho;
    }

    distEPhi=std::fabs(pPhi-fSPhi-fDPhi)*rho;

    // Find new minimum
    //
    if (distSPhi<distEPhi)
    {
      if (distSPhi<distMin)
      {
        distMin=distSPhi;
        side=kNSPhi;
      }
    }
    else
    {
      if (distEPhi<distMin)
      {
        distMin=distEPhi;
        side=kNEPhi;
      }
    }
  }

  //
  // Distance to theta planes
  //

  if (!fFullThetaSphere && radius)
  {
    pTheta=std::atan2(rho,p.z());
    distSTheta=std::fabs(pTheta-fSTheta)*radius;
    distETheta=std::fabs(pTheta-fSTheta-fDTheta)*radius;

    // Find new minimum
    //
    if (distSTheta<distETheta)
    {
      if (distSTheta<distMin)
      {
        distMin = distSTheta ;
        side = kNSTheta ;
      }
    }
    else
    {
      if (distETheta<distMin)
      {
        distMin = distETheta ;
        side = kNETheta ;
      }
    }
  }

  switch (side)
  {
    case kNRMin:      // Inner radius
      norm=G4ThreeVector(-p.x()/radius,-p.y()/radius,-p.z()/radius);
      break;
    case kNRMax:      // Outer radius
      norm=G4ThreeVector(p.x()/radius,p.y()/radius,p.z()/radius);
      break;
    case kNSPhi:
      norm=G4ThreeVector(sinSPhi,-cosSPhi,0);
      break;
    case kNEPhi:
      norm=G4ThreeVector(-sinEPhi,cosEPhi,0);
      break;
    case kNSTheta:
      norm=G4ThreeVector(-cosSTheta*std::cos(pPhi),
                         -cosSTheta*std::sin(pPhi),
                          sinSTheta            );
      break;
    case kNETheta:
      norm=G4ThreeVector( cosETheta*std::cos(pPhi),
                          cosETheta*std::sin(pPhi),
                         -sinETheta              );
      break;
    default:          // Should never reach this case ...
      DumpInfo();
      G4Exception("G4Sphere::ApproxSurfaceNormal()",
                  "GeomSolids1002", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break;    
  }

  return norm;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// -> If point is outside outer radius, compute intersection with rmax
//        - if no intersection return
//        - if  valid phi,theta return intersection Dist
//
// -> If shell, compute intersection with inner radius, taking largest +ve root
//        - if valid phi,theta, save intersection
//
// -> If phi segmented, compute intersection with phi half planes
//        - if valid intersection(r,theta), return smallest intersection of
//          inner shell & phi intersection
//
// -> If theta segmented, compute intersection with theta cones
//        - if valid intersection(r,phi), return smallest intersection of
//          inner shell & theta intersection
//
//
// NOTE:
// - `if valid' (above) implies tolerant checking of intersection points
//
// OPT:
// Move tolIO/ORmin/RMax2 precalcs to where they are needed -
// not required for most cases.
// Avoid atan2 for non theta cut G4Sphere.

G4double G4Sphere::DistanceToIn( const G4ThreeVector& p,
                                 const G4ThreeVector& v  ) const
{
  G4double snxt = kInfinity ;      // snxt = default return value
  G4double rho2, rad2, pDotV2d, pDotV3d, pTheta ;
  G4double tolSTheta=0., tolETheta=0. ;
  const G4double dRmax = 100.*fRmax;

  const G4double halfRmaxTolerance = fRmaxTolerance*0.5;
  const G4double halfRminTolerance = fRminTolerance*0.5;
  const G4double tolORMin2 = (fRmin>halfRminTolerance)
               ? (fRmin-halfRminTolerance)*(fRmin-halfRminTolerance) : 0;
  const G4double tolIRMin2 =
               (fRmin+halfRminTolerance)*(fRmin+halfRminTolerance);
  const G4double tolORMax2 =
               (fRmax+halfRmaxTolerance)*(fRmax+halfRmaxTolerance);
  const G4double tolIRMax2 =
               (fRmax-halfRmaxTolerance)*(fRmax-halfRmaxTolerance);

  // Intersection point
  //
  G4double xi, yi, zi, rhoi, rhoi2, radi2, iTheta ;

  // Phi intersection
  //
  G4double Comp ; 

  // Phi precalcs
  //
  G4double Dist, cosPsi ;

  // Theta precalcs
  //
  G4double dist2STheta, dist2ETheta ;
  G4double t1, t2, b, c, d2, d, sd = kInfinity ;

  // General Precalcs
  //
  rho2 = p.x()*p.x() + p.y()*p.y() ;
  rad2 = rho2 + p.z()*p.z() ;
  pTheta = std::atan2(std::sqrt(rho2),p.z()) ;

  pDotV2d = p.x()*v.x() + p.y()*v.y() ;
  pDotV3d = pDotV2d + p.z()*v.z() ;

  // Theta precalcs
  //
  if (!fFullThetaSphere)
  {
    tolSTheta = fSTheta - halfAngTolerance ;
    tolETheta = eTheta + halfAngTolerance ;

    // Special case rad2 = 0 comparing with direction
    //
    if ((rad2!=0.0) || (fRmin!=0.0))
    {
      // Keep going for computation of distance...
    }
    else  // Positioned on the sphere's origin
    {
      G4double vTheta = std::atan2(std::sqrt(v.x()*v.x()+v.y()*v.y()),v.z()) ;
      if ( (vTheta < tolSTheta) || (vTheta > tolETheta) )
      {
        return snxt ; // kInfinity
      }
      return snxt = 0.0 ;
    }
  }

  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer G4Sphere heading through solid (-> 0)
  // - No intersect -> no intersection with G4Sphere
  //
  // Shell eqn: x^2+y^2+z^2=RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2sd(pxvx+pyvy+pzvz)+sd^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2sd(pDotV3d)       +sd^2                =R^2
  //
  // => sd=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))

  c = rad2 - fRmax*fRmax ;

  if (c > fRmaxTolerance*fRmax)
  {
    // If outside tolerant boundary of outer G4Sphere
    // [should be std::sqrt(rad2)-fRmax > halfRmaxTolerance]

    d2 = pDotV3d*pDotV3d - c ;

    if ( d2 >= 0 )
    {
      sd = -pDotV3d - std::sqrt(d2) ;

      if (sd >= 0 )
      {
        if ( sd>dRmax ) // Avoid rounding errors due to precision issues seen on
        {               // 64 bits systems. Split long distances and recompute
          G4double fTerm = sd-std::fmod(sd,dRmax);
          sd = fTerm + DistanceToIn(p+fTerm*v,v);
        } 
        xi   = p.x() + sd*v.x() ;
        yi   = p.y() + sd*v.y() ;
        rhoi = std::sqrt(xi*xi + yi*yi) ;

        if (!fFullPhiSphere && rhoi)    // Check phi intersection
        {
          cosPsi = (xi*cosCPhi + yi*sinCPhi)/rhoi ;

          if (cosPsi >= cosHDPhiOT)
          {
            if (!fFullThetaSphere)   // Check theta intersection
            {
              zi = p.z() + sd*v.z() ;

              // rhoi & zi can never both be 0
              // (=>intersect at origin =>fRmax=0)
              //
              iTheta = std::atan2(rhoi,zi) ;
              if ( (iTheta >= tolSTheta) && (iTheta <= tolETheta) )
              {
                return snxt = sd ;
              }
            }
            else
            {
              return snxt=sd;
            }
          }
        }
        else
        {
          if (!fFullThetaSphere)    // Check theta intersection
          {
            zi = p.z() + sd*v.z() ;

            // rhoi & zi can never both be 0
            // (=>intersect at origin => fRmax=0 !)
            //
            iTheta = std::atan2(rhoi,zi) ;
            if ( (iTheta >= tolSTheta) && (iTheta <= tolETheta) )
            {
              return snxt=sd;
            }
          }
          else
          {
            return snxt = sd;
          }
        }          
      }
    }
    else    // No intersection with G4Sphere
    {
      return snxt=kInfinity;
    }
  }
  else
  {
    // Inside outer radius
    // check not inside, and heading through G4Sphere (-> 0 to in)

    d2 = pDotV3d*pDotV3d - c ;

    if ( (rad2 > tolIRMax2)
      && ( (d2 >= fRmaxTolerance*fRmax) && (pDotV3d < 0) ) )
    {
      if (!fFullPhiSphere)
      {
        // Use inner phi tolerant boundary -> if on tolerant
        // phi boundaries, phi intersect code handles leaving/entering checks

        cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/std::sqrt(rho2) ;

        if (cosPsi>=cosHDPhiIT)
        { 
          // inside radii, delta r -ve, inside phi

          if ( !fFullThetaSphere )
          {
            if ( (pTheta >= tolSTheta + kAngTolerance)
              && (pTheta <= tolETheta - kAngTolerance) )
            {
              return snxt=0;
            }
          }
          else    // strictly inside Theta in both cases
          {
            return snxt=0;
          }
        }
      }
      else
      {
        if ( !fFullThetaSphere )
        {
          if ( (pTheta >= tolSTheta + kAngTolerance)
            && (pTheta <= tolETheta - kAngTolerance) )
          {
            return snxt=0;
          }
        }
        else   // strictly inside Theta in both cases
        {
          return snxt=0;
        }
      }
    }
  }

  // Inner spherical shell intersection
  // - Always farthest root, because would have passed through outer
  //   surface first.
  // - Tolerant check if travelling through solid

  if (fRmin)
  {
    c  = rad2 - fRmin*fRmin ;
    d2 = pDotV3d*pDotV3d - c ;

    // Within tolerance inner radius of inner G4Sphere
    // Check for immediate entry/already inside and travelling outwards

    if ( (c > -halfRminTolerance) && (rad2 < tolIRMin2)
      && ( (d2 < fRmin*kCarTolerance) || (pDotV3d >= 0) ) )
    {
      if ( !fFullPhiSphere )
      {
        // Use inner phi tolerant boundary -> if on tolerant
        // phi boundaries, phi intersect code handles leaving/entering checks

        cosPsi = (p.x()*cosCPhi+p.y()*sinCPhi)/std::sqrt(rho2) ;
        if (cosPsi >= cosHDPhiIT)
        { 
          // inside radii, delta r -ve, inside phi
          //
          if ( !fFullThetaSphere )
          {
            if ( (pTheta >= tolSTheta + kAngTolerance)
              && (pTheta <= tolETheta - kAngTolerance) )
            {
              return snxt=0;
            }
          }
          else
          {
            return snxt = 0 ;
          }
        }
      }
      else
      {
        if ( !fFullThetaSphere )
        {
          if ( (pTheta >= tolSTheta + kAngTolerance)
            && (pTheta <= tolETheta - kAngTolerance) )
          {
            return snxt = 0 ;
          }
        }
        else
        {
          return snxt=0;
        }
      }
    }
    else   // Not special tolerant case
    {
      if (d2 >= 0)
      {
        sd = -pDotV3d + std::sqrt(d2) ;
        if ( sd >= halfRminTolerance )  // It was >= 0 ??
        {
          xi   = p.x() + sd*v.x() ;
          yi   = p.y() + sd*v.y() ;
          rhoi = std::sqrt(xi*xi+yi*yi) ;

          if ( !fFullPhiSphere && rhoi )   // Check phi intersection
          {
            cosPsi = (xi*cosCPhi + yi*sinCPhi)/rhoi ;

            if (cosPsi >= cosHDPhiOT)
            {
              if ( !fFullThetaSphere )  // Check theta intersection
              {
                zi = p.z() + sd*v.z() ;

                // rhoi & zi can never both be 0
                // (=>intersect at origin =>fRmax=0)
                //
                iTheta = std::atan2(rhoi,zi) ;
                if ( (iTheta >= tolSTheta) && (iTheta<=tolETheta) )
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt=sd;
              }
            }
          }
          else
          {
            if ( !fFullThetaSphere )   // Check theta intersection
            {
              zi = p.z() + sd*v.z() ;

              // rhoi & zi can never both be 0
              // (=>intersect at origin => fRmax=0 !)
              //
              iTheta = std::atan2(rhoi,zi) ;
              if ( (iTheta >= tolSTheta) && (iTheta <= tolETheta) )
              {
                snxt = sd;
              }
            }
            else
            {
              snxt = sd;
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
  //
  if ( !fFullPhiSphere )
  {
    // First phi surface ('S'tarting phi)
    // Comp = Component in outwards normal dirn
    //
    Comp = v.x()*sinSPhi - v.y()*cosSPhi ;
                    
    if ( Comp < 0 )
    {
      Dist = p.y()*cosSPhi - p.x()*sinSPhi ;

      if (Dist < halfCarTolerance)
      {
        sd = Dist/Comp ;

        if (sd < snxt)
        {
          if ( sd > 0 )
          {
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;
          }
          else
          {
            sd    = 0     ;
            xi    = p.x() ;
            yi    = p.y() ;
            zi    = p.z() ;
            rhoi2 = rho2  ;
            radi2 = rad2  ;
          }
          if ( (radi2 <= tolORMax2)
            && (radi2 >= tolORMin2)
            && ((yi*cosCPhi-xi*sinCPhi) <= 0) )
          {
            // Check theta intersection
            // rhoi & zi can never both be 0
            // (=>intersect at origin =>fRmax=0)
            //
            if ( !fFullThetaSphere )
            {
              iTheta = std::atan2(std::sqrt(rhoi2),zi) ;
              if ( (iTheta >= tolSTheta) && (iTheta <= tolETheta) )
              {
                // r and theta intersections good
                // - check intersecting with correct half-plane

                if ((yi*cosCPhi-xi*sinCPhi) <= 0)
                {
                  snxt = sd;
                }
              }
            }
            else
            {
              snxt = sd;
            }
          }
        }
      }
    }

    // Second phi surface ('E'nding phi)
    // Component in outwards normal dirn

    Comp = -( v.x()*sinEPhi-v.y()*cosEPhi ) ;
        
    if (Comp < 0)
    {
      Dist = -(p.y()*cosEPhi-p.x()*sinEPhi) ;
      if ( Dist < halfCarTolerance )
      {
        sd = Dist/Comp ;

        if ( sd < snxt )
        {
          if (sd > 0)
          {
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;
          }
          else
          {
            sd    = 0     ;
            xi    = p.x() ;
            yi    = p.y() ;
            zi    = p.z() ;
            rhoi2 = rho2  ;
            radi2 = rad2  ;
          }
          if ( (radi2 <= tolORMax2)
            && (radi2 >= tolORMin2)
            && ((yi*cosCPhi-xi*sinCPhi) >= 0) )
          {
            // Check theta intersection
            // rhoi & zi can never both be 0
            // (=>intersect at origin =>fRmax=0)
            //
            if ( !fFullThetaSphere )
            {
              iTheta = std::atan2(std::sqrt(rhoi2),zi) ;
              if ( (iTheta >= tolSTheta) && (iTheta <= tolETheta) )
              {
                // r and theta intersections good
                // - check intersecting with correct half-plane

                if ((yi*cosCPhi-xi*sinCPhi) >= 0)
                {
                  snxt = sd;
                }
              }
            }
            else
            {
              snxt = sd;
            }
          }
        }
      }
    }
  }

  // Theta segment intersection

  if ( !fFullThetaSphere )
  {

    // Intersection with theta surfaces
    // Known failure cases:
    // o  Inside tolerance of stheta surface, skim
    //    ~parallel to cone and Hit & enter etheta surface [& visa versa]
    //
    //    To solve: Check 2nd root of etheta surface in addition to stheta
    //
    // o  start/end theta is exactly pi/2 
    // Intersections with cones
    //
    // Cone equation: x^2+y^2=z^2tan^2(t)
    //
    // => (px+svx)^2+(py+svy)^2=(pz+svz)^2tan^2(t)
    //
    // => (px^2+py^2-pz^2tan^2(t))+2sd(pxvx+pyvy-pzvztan^2(t))
    //       + sd^2(vx^2+vy^2-vz^2tan^2(t)) = 0
    //
    // => sd^2(1-vz^2(1+tan^2(t))+2sd(pdotv2d-pzvztan^2(t))
    //       + (rho2-pz^2tan^2(t)) = 0

    if (fSTheta)
    {
      dist2STheta = rho2 - p.z()*p.z()*tanSTheta2 ;
    }
    else
    {
      dist2STheta = kInfinity ;
    }
    if ( eTheta < pi )
    {
      dist2ETheta=rho2-p.z()*p.z()*tanETheta2;
    }
    else
    {
      dist2ETheta=kInfinity;
    }      
    if ( pTheta < tolSTheta )
    {
      // Inside (theta<stheta-tol) stheta cone
      // First root of stheta cone, second if first root -ve

      t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
      t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;
      if (t1)
      {   
        b  = t2/t1 ;
        c  = dist2STheta/t1 ;
        d2 = b*b - c ;

        if ( d2 >= 0 )
        {
          d  = std::sqrt(d2) ;
          sd = -b - d ;    // First root
          zi = p.z() + sd*v.z();

          if ( (sd < 0) || (zi*(fSTheta - halfpi) > 0) )
          {
            sd = -b+d;    // Second root
          }
          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd*v.x();
            yi    = p.y() + sd*v.y();
            zi    = p.z() + sd*v.z();
            rhoi2 = xi*xi + yi*yi;
            radi2 = rhoi2 + zi*zi;
            if ( (radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && (zi*(fSTheta - halfpi) <= 0) )
            {
              if ( !fFullPhiSphere && rhoi2 )  // Check phi intersection
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }

      // Possible intersection with ETheta cone. 
      // Second >= 0 root should be considered
        
      if ( eTheta < pi )
      {
        t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
        t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;
        if (t1)
        { 
          b  = t2/t1 ;
          c  = dist2ETheta/t1 ;
          d2 = b*b - c ;

          if (d2 >= 0)
          {
            d  = std::sqrt(d2) ;
            sd = -b + d ;    // Second root

            if ( (sd >= 0) && (sd < snxt) )
            {
              xi    = p.x() + sd*v.x() ;
              yi    = p.y() + sd*v.y() ;
              zi    = p.z() + sd*v.z() ;
              rhoi2 = xi*xi + yi*yi   ;
              radi2 = rhoi2 + zi*zi   ;

              if ( (radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi*(eTheta - halfpi) <= 0) )
              {
                if (!fFullPhiSphere && rhoi2)   // Check phi intersection
                {
                  cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                  if (cosPsi >= cosHDPhiOT)
                  {
                    snxt = sd;
                  }
                }
                else
                {
                  snxt = sd;
                }
              }
            }
          }
        }
      }
    }  
    else if ( pTheta > tolETheta ) 
    { 
      // dist2ETheta<-kRadTolerance*0.5 && dist2STheta>0)
      // Inside (theta > etheta+tol) e-theta cone
      // First root of etheta cone, second if first root 'imaginary'

      t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
      t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;
      if (t1)
      {  
        b  = t2/t1 ;
        c  = dist2ETheta/t1 ;
        d2 = b*b - c ;

        if (d2 >= 0)
        {
          d  = std::sqrt(d2) ;
          sd = -b - d ;    // First root
          zi = p.z() + sd*v.z();

          if ( (sd < 0) || (zi*(eTheta - halfpi) > 0) )
          {
            sd = -b + d ;           // second root
          }
          if ( (sd >= 0) && (sd < snxt) )
          {
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;

            if ( (radi2 <= tolORMax2)
              && (radi2 >= tolORMin2) 
              && (zi*(eTheta - halfpi) <= 0) )
            {
              if (!fFullPhiSphere && rhoi2)  // Check phi intersection
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }

      // Possible intersection with STheta cone. 
      // Second >= 0 root should be considered
        
      if ( fSTheta )
      {
        t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
        t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;
        if (t1)
        { 
          b  = t2/t1 ;
          c  = dist2STheta/t1 ;
          d2 = b*b - c ;

          if (d2 >= 0)
          {
            d  = std::sqrt(d2) ;
            sd = -b + d ;    // Second root

            if ( (sd >= 0) && (sd < snxt) )
            {
              xi    = p.x() + sd*v.x() ;
              yi    = p.y() + sd*v.y() ;
              zi    = p.z() + sd*v.z() ;
              rhoi2 = xi*xi + yi*yi   ;
              radi2 = rhoi2 + zi*zi   ;

              if ( (radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi*(fSTheta - halfpi) <= 0) )
              {
                if (!fFullPhiSphere && rhoi2)   // Check phi intersection
                {
                  cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                  if (cosPsi >= cosHDPhiOT)
                  {
                    snxt = sd;
                  }
                }
                else
                {
                  snxt = sd;
                }
              }
            }
          }
        }
      }  
    }     
    else if ( (pTheta < tolSTheta + kAngTolerance)
           && (fSTheta > halfAngTolerance) )
    {
      // In tolerance of stheta
      // If entering through solid [r,phi] => 0 to in
      // else try 2nd root

      t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;
      if ( (t2>=0 && tolIRMin2<rad2 && rad2<tolIRMax2 && fSTheta<halfpi)
        || (t2<0  && tolIRMin2<rad2 && rad2<tolIRMax2 && fSTheta>halfpi)
        || (v.z()<0 && tolIRMin2<rad2 && rad2<tolIRMax2 && fSTheta==halfpi) )
      {
        if (!fFullPhiSphere && rho2)  // Check phi intersection
        {
          cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/std::sqrt(rho2) ;
          if (cosPsi >= cosHDPhiIT)
          {
            return 0 ;
          }
        }
        else
        {
          return 0 ;
        }
      }

      // Not entering immediately/travelling through

      t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
      if (t1)
      { 
        b  = t2/t1 ;
        c  = dist2STheta/t1 ;
        d2 = b*b - c ;

        if (d2 >= 0)
        {
          d  = std::sqrt(d2) ;
          sd = -b + d ;
          if ( (sd >= halfCarTolerance) && (sd < snxt) && (fSTheta < halfpi) )
          {   // ^^^^^^^^^^^^^^^^^^^^^  shouldn't it be >=0 instead ?
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;

            if ( (radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && (zi*(fSTheta - halfpi) <= 0) )
            {
              if ( !fFullPhiSphere && rhoi2 )    // Check phi intersection
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                if ( cosPsi >= cosHDPhiOT )
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }   
    else if ((pTheta > tolETheta-kAngTolerance) && (eTheta < pi-kAngTolerance))
    {

      // In tolerance of etheta
      // If entering through solid [r,phi] => 0 to in
      // else try 2nd root

      t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;

      if (   ((t2<0) && (eTheta < halfpi)
          && (tolIRMin2 < rad2) && (rad2 < tolIRMax2))
        ||   ((t2>=0) && (eTheta > halfpi)
          && (tolIRMin2 < rad2) && (rad2 < tolIRMax2))
        ||   ((v.z()>0) && (eTheta == halfpi)
          && (tolIRMin2 < rad2) && (rad2 < tolIRMax2))  )
      {
        if (!fFullPhiSphere && rho2)   // Check phi intersection
        {
          cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/std::sqrt(rho2) ;
          if (cosPsi >= cosHDPhiIT)
          {
            return 0 ;
          }
        }
        else
        {
          return 0 ;
        }
      }

      // Not entering immediately/travelling through

      t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
      if (t1)
      { 
        b  = t2/t1 ;
        c  = dist2ETheta/t1 ;
        d2 = b*b - c ;

        if (d2 >= 0)
        {
          d  = std::sqrt(d2) ;
          sd = -b + d ;
        
          if ( (sd >= halfCarTolerance)
            && (sd < snxt) && (eTheta > halfpi) )
          {
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;

            if ( (radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && (zi*(eTheta - halfpi) <= 0) )
            {
              if (!fFullPhiSphere && rhoi2)   // Check phi intersection
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        } 
      }   
    }  
    else
    {
      // stheta+tol<theta<etheta-tol
      // For BOTH stheta & etheta check 2nd root for validity [r,phi]

      t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
      t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;
      if (t1)
      { 
        b  = t2/t1;
        c  = dist2STheta/t1 ;
        d2 = b*b - c ;

        if (d2 >= 0)
        {
          d  = std::sqrt(d2) ;
          sd = -b + d ;    // second root

          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;

            if ( (radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && (zi*(fSTheta - halfpi) <= 0) )
            {
              if (!fFullPhiSphere && rhoi2)   // Check phi intersection
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }        
      t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
      t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;
      if (t1)
      {   
        b  = t2/t1 ;
        c  = dist2ETheta/t1 ;
        d2 = b*b - c ;

        if (d2 >= 0)
        {
          d  = std::sqrt(d2) ;
          sd = -b + d;    // second root

          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd*v.x() ;
            yi    = p.y() + sd*v.y() ;
            zi    = p.z() + sd*v.z() ;
            rhoi2 = xi*xi + yi*yi   ;
            radi2 = rhoi2 + zi*zi   ;

            if ( (radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && (zi*(eTheta - halfpi) <= 0) )
            {
              if (!fFullPhiSphere && rhoi2)   // Check phi intersection
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/std::sqrt(rhoi2) ;
                if ( cosPsi >= cosHDPhiOT )
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }  
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to radial planes
// - Only to phi planes if outside phi extent
// - Only to theta planes if outside theta extent
// - Return 0 if point inside

G4double G4Sphere::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe=0.0,safeRMin,safeRMax,safePhi,safeTheta;
  G4double rho2,rds,rho;
  G4double cosPsi;
  G4double pTheta,dTheta1,dTheta2;
  rho2=p.x()*p.x()+p.y()*p.y();
  rds=std::sqrt(rho2+p.z()*p.z());
  rho=std::sqrt(rho2);

  //
  // Distance to r shells
  //    
  if (fRmin)
  {
    safeRMin=fRmin-rds;
    safeRMax=rds-fRmax;
    if (safeRMin>safeRMax)
    {
      safe=safeRMin;
    }
    else
    {
      safe=safeRMax;
    }
  }
  else
  {
    safe=rds-fRmax;
  }

  //
  // Distance to phi extent
  //
  if (!fFullPhiSphere && rho)
  {
    // Psi=angle from central phi to point
    //
    cosPsi=(p.x()*cosCPhi+p.y()*sinCPhi)/rho;
    if (cosPsi<std::cos(hDPhi))
    {
      // Point lies outside phi range
      //
      if ((p.y()*cosCPhi-p.x()*sinCPhi)<=0)
      {
        safePhi=std::fabs(p.x()*sinSPhi-p.y()*cosSPhi);
      }
      else
      {
        safePhi=std::fabs(p.x()*sinEPhi-p.y()*cosEPhi);
      }
      if (safePhi>safe)  { safe=safePhi; }
    }
  }
  //
  // Distance to Theta extent
  //    
  if ((rds!=0.0) && (!fFullThetaSphere))
  {
    pTheta=std::acos(p.z()/rds);
    if (pTheta<0)  { pTheta+=pi; }
    dTheta1=fSTheta-pTheta;
    dTheta2=pTheta-eTheta;
    if (dTheta1>dTheta2)
    {
      if (dTheta1>=0)             // WHY ???????????
      {
        safeTheta=rds*std::sin(dTheta1);
        if (safe<=safeTheta)
        {
          safe=safeTheta;
        }
      }
    }
    else
    {
      if (dTheta2>=0)
      {
        safeTheta=rds*std::sin(dTheta2);
        if (safe<=safeTheta)
        {
          safe=safeTheta;
        }
      }
    }
  }

  if (safe<0)  { safe=0; }
  return safe;
}

/////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from 'inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4Sphere::DistanceToOut( const G4ThreeVector& p,
                                  const G4ThreeVector& v,
                                  const G4bool calcNorm,
                                        G4bool *validNorm,
                                        G4ThreeVector *n   ) const
{
  G4double snxt = kInfinity;     // snxt is default return value
  G4double sphi= kInfinity,stheta= kInfinity;
  ESide side=kNull,sidephi=kNull,sidetheta=kNull;  

  const G4double halfRmaxTolerance = fRmaxTolerance*0.5;
  const G4double halfRminTolerance = fRminTolerance*0.5;
  const G4double Rmax_plus  = fRmax + halfRmaxTolerance;
  const G4double Rmin_minus = (fRmin) ? fRmin-halfRminTolerance : 0;
  G4double t1,t2;
  G4double b,c,d;

  // Variables for phi intersection:

  G4double pDistS,compS,pDistE,compE,sphi2,vphi;
    
  G4double rho2,rad2,pDotV2d,pDotV3d;

  G4double xi,yi,zi;      // Intersection point

  // Theta precals
  //
  G4double rhoSecTheta;
  G4double dist2STheta, dist2ETheta, distTheta;
  G4double d2,sd;

  // General Precalcs
  //
  rho2 = p.x()*p.x()+p.y()*p.y();
  rad2 = rho2+p.z()*p.z();

  pDotV2d = p.x()*v.x()+p.y()*v.y();
  pDotV3d = pDotV2d+p.z()*v.z();

  // Radial Intersections from G4Sphere::DistanceToIn
  //
  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer G4Sphere heading through solid (-> 0)
  // - No intersect -> no intersection with G4Sphere
  //
  // Shell eqn: x^2+y^2+z^2=RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2sd(pxvx+pyvy+pzvz)+sd^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2sd(pDotV3d)       +sd^2                =R^2
  //
  // => sd=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))

  if( (rad2 <= Rmax_plus*Rmax_plus) && (rad2 >= Rmin_minus*Rmin_minus) )
  {
    c = rad2 - fRmax*fRmax;

    if (c < fRmaxTolerance*fRmax) 
    {
      // Within tolerant Outer radius 
      // 
      // The test is
      //     rad  - fRmax < 0.5*kRadTolerance
      // =>  rad  < fRmax + 0.5*kRadTol
      // =>  rad2 < (fRmax + 0.5*kRadTol)^2
      // =>  rad2 < fRmax^2 + 2.*0.5*fRmax*kRadTol + 0.25*kRadTol*kRadTol
      // =>  rad2 - fRmax^2    <~    fRmax*kRadTol 

      d2 = pDotV3d*pDotV3d - c;

      if( (c >- fRmaxTolerance*fRmax)       // on tolerant surface
       && ((pDotV3d >=0) || (d2 < 0)) )     // leaving outside from Rmax 
                                            // not re-entering
      {
        if(calcNorm)
        {
          *validNorm = true ;
          *n         = G4ThreeVector(p.x()/fRmax,p.y()/fRmax,p.z()/fRmax) ;
        }
        return snxt = 0;
      }
      else 
      {
        snxt = -pDotV3d+std::sqrt(d2);    // second root since inside Rmax
        side =  kRMax ; 
      }
    }

    // Inner spherical shell intersection:
    // Always first >=0 root, because would have passed
    // from outside of Rmin surface .

    if (fRmin)
    {
      c  = rad2 - fRmin*fRmin;
      d2 = pDotV3d*pDotV3d - c;

      if (c >- fRminTolerance*fRmin) // 2.0 * (0.5*kRadTolerance) * fRmin
      {
        if ( (c < fRminTolerance*fRmin)              // leaving from Rmin
          && (d2 >= fRminTolerance*fRmin) && (pDotV3d < 0) )
        {
          if(calcNorm)  { *validNorm = false; }  // Rmin surface is concave
          return snxt = 0 ;
        }
        else
        {  
          if ( d2 >= 0. )
          {
            sd = -pDotV3d-std::sqrt(d2);

            if ( sd >= 0. )     // Always intersect Rmin first
            {
              snxt = sd ;
              side = kRMin ;
            }
          }
        }
      }
    }
  }

  // Theta segment intersection

  if ( !fFullThetaSphere )
  {
    // Intersection with theta surfaces
    //
    // Known failure cases:
    // o  Inside tolerance of stheta surface, skim
    //    ~parallel to cone and Hit & enter etheta surface [& visa versa]
    //
    //    To solve: Check 2nd root of etheta surface in addition to stheta
    //
    // o  start/end theta is exactly pi/2 
    //
    // Intersections with cones
    //
    // Cone equation: x^2+y^2=z^2tan^2(t)
    //
    // => (px+svx)^2+(py+svy)^2=(pz+svz)^2tan^2(t)
    //
    // => (px^2+py^2-pz^2tan^2(t))+2sd(pxvx+pyvy-pzvztan^2(t))
    //       + sd^2(vx^2+vy^2-vz^2tan^2(t)) = 0
    //
    // => sd^2(1-vz^2(1+tan^2(t))+2sd(pdotv2d-pzvztan^2(t))
    //       + (rho2-pz^2tan^2(t)) = 0
    //
  
    if(fSTheta) // intersection with first cons
    {
      if( std::fabs(tanSTheta) > 5./kAngTolerance ) // kons is plane z=0
      {
        if( v.z() > 0. ) 
        {
          if ( std::fabs( p.z() ) <= halfRmaxTolerance )
          {
            if(calcNorm)
            {
              *validNorm = true;
              *n = G4ThreeVector(0.,0.,1.);
            }
            return snxt = 0 ;
          }  
          stheta    = -p.z()/v.z();
          sidetheta = kSTheta;
        }
      }
      else // kons is not plane 
      {
        t1          = 1-v.z()*v.z()*(1+tanSTheta2);
        t2          = pDotV2d-p.z()*v.z()*tanSTheta2;  // ~vDotN if p on cons
        dist2STheta = rho2-p.z()*p.z()*tanSTheta2;     // t3

        distTheta = std::sqrt(rho2)-p.z()*tanSTheta;

        if( std::fabs(t1) < halfAngTolerance ) // 1st order equation,
        {                                      // v parallel to kons
          if( v.z() > 0. )
          {
            if(std::fabs(distTheta) < halfRmaxTolerance) // p on surface
            {
              if( (fSTheta < halfpi) && (p.z() > 0.) )
              {
                if( calcNorm )  { *validNorm = false; }
                return snxt = 0.;
              }
              else if( (fSTheta > halfpi) && (p.z() <= 0) )
              {
                if( calcNorm ) 
                {
                  *validNorm = true;
                  if (rho2)
                  {
                    rhoSecTheta = std::sqrt(rho2*(1+tanSTheta2));
                   
                    *n = G4ThreeVector( p.x()/rhoSecTheta,   
                                        p.y()/rhoSecTheta,
                                        std::sin(fSTheta)  );
                  }
                  else *n = G4ThreeVector(0.,0.,1.);
                }
                return snxt = 0.;               
              }
            }
            stheta    = -0.5*dist2STheta/t2;
            sidetheta = kSTheta;
          }  
        }      // 2nd order equation, 1st root of fSTheta cone,
        else   // 2nd if 1st root -ve
        {
          if( std::fabs(distTheta) < halfRmaxTolerance )
          {
            if( (fSTheta > halfpi) && (t2 >= 0.) ) // leave
            {
              if( calcNorm ) 
              {
                *validNorm = true;
                if (rho2)
                {
                  rhoSecTheta = std::sqrt(rho2*(1+tanSTheta2));
                   
                  *n = G4ThreeVector( p.x()/rhoSecTheta,   
                                      p.y()/rhoSecTheta,
                                      std::sin(fSTheta)  );
                }
                else  { *n = G4ThreeVector(0.,0.,1.); }
              }
              return snxt = 0.;
            }
            else if( (fSTheta < halfpi) && (t2 < 0.) && (p.z() >=0.) ) // leave
            {
              if( calcNorm )  { *validNorm = false; }
              return snxt = 0.;
            }                               
          }
          b  = t2/t1;
          c  = dist2STheta/t1;
          d2 = b*b - c ;

          if ( d2 >= 0. )
          {
            d = std::sqrt(d2);

            if( fSTheta > halfpi )
            {
              sd = -b - d;         // First root

              if ( ((std::fabs(s) < halfRmaxTolerance) && (t2 < 0.))
               ||  (sd < 0.)  || ( (sd > 0.) && (p.z() + sd*v.z() > 0.) )     ) 
              {
                sd = -b + d ; // 2nd root
              }
              if( (sd > halfRmaxTolerance) && (p.z() + sd*v.z() <= 0.) )  
              {
                stheta    = sd;
                sidetheta = kSTheta;
              }
            }
            else // sTheta < pi/2, concave surface, no normal
            {
              sd = -b - d;         // First root

              if ( ( (std::fabs(sd) < halfRmaxTolerance) && (t2 >= 0.) )
                || (sd < 0.) || ( (sd > 0.) && (p.z() + sd*v.z() < 0.) )   )
              {
                sd = -b + d ; // 2nd root
              }
              if( (sd > halfRmaxTolerance) && (p.z() + sd*v.z() >= 0.) )  
              {
                stheta    = sd;
                sidetheta = kSTheta;
              }            
            }
          }
        }
      }
    }
    if (eTheta < pi) // intersection with second cons
    {
      if( std::fabs(tanETheta) > 5./kAngTolerance ) // kons is plane z=0
      {
        if( v.z() < 0. ) 
        {
          if ( std::fabs( p.z() ) <= halfRmaxTolerance )
          {
            if(calcNorm)
            {
              *validNorm = true;
              *n = G4ThreeVector(0.,0.,-1.);
            }
            return snxt = 0 ;
          }  
          sd = -p.z()/v.z();

          if( sd < stheta )
          {
            stheta    = sd;
            sidetheta = kETheta;
          }
        }
      }
      else // kons is not plane 
      {
        t1          = 1-v.z()*v.z()*(1+tanETheta2);
        t2          = pDotV2d-p.z()*v.z()*tanETheta2;  // ~vDotN if p on cons
        dist2ETheta = rho2-p.z()*p.z()*tanETheta2;     // t3

        distTheta = std::sqrt(rho2)-p.z()*tanETheta;

        if( std::fabs(t1) < halfAngTolerance ) // 1st order equation,
        {                                      // v parallel to kons
          if( v.z() < 0. )
          {
            if(std::fabs(distTheta) < halfRmaxTolerance) // p on surface
            {
              if( (eTheta > halfpi) && (p.z() < 0.) )
              {
                if( calcNorm )  { *validNorm = false; }
                return snxt = 0.;
              }
              else if ( (eTheta < halfpi) && (p.z() >= 0) )
              {
                if( calcNorm ) 
                {
                  *validNorm = true;
                  if (rho2)
                  {
                    rhoSecTheta = std::sqrt(rho2*(1+tanETheta2));
                    *n = G4ThreeVector( p.x()/rhoSecTheta,   
                                        p.y()/rhoSecTheta,
                                        -sinETheta  );
                  }
                  else  { *n = G4ThreeVector(0.,0.,-1.); }
                }
                return snxt = 0.;               
              }
            }
            sd = -0.5*dist2ETheta/t2;

            if( sd < stheta )
            {
              stheta    = sd;
              sidetheta = kETheta;
            }
          }  
        }      // 2nd order equation, 1st root of fSTheta cone
        else   // 2nd if 1st root -ve
        {
          if ( std::fabs(distTheta) < halfRmaxTolerance )
          {
            if( (eTheta < halfpi) && (t2 >= 0.) ) // leave
            {
              if( calcNorm ) 
              {
                *validNorm = true;
                if (rho2)
                {
                    rhoSecTheta = std::sqrt(rho2*(1+tanETheta2));
                    *n = G4ThreeVector( p.x()/rhoSecTheta,   
                                        p.y()/rhoSecTheta,
                                        -sinETheta  );
                }
                else *n = G4ThreeVector(0.,0.,-1.);
              }                           
              return snxt = 0.;
            }
            else if ( (eTheta > halfpi)
                   && (t2 < 0.) && (p.z() <=0.) ) // leave
            {
              if( calcNorm )  { *validNorm = false; }
              return snxt = 0.;
            }                               
          }
          b  = t2/t1;
          c  = dist2ETheta/t1;
          d2 = b*b - c ;
          if ( (d2 <halfRmaxTolerance) && (d2 > -halfRmaxTolerance) )
          {
            d2 = 0.;
          }
          if ( d2 >= 0. )
          {
            d = std::sqrt(d2);

            if( eTheta < halfpi )
            {
              sd = -b - d;         // First root

              if( ((std::fabs(sd) < halfRmaxTolerance) && (t2 < 0.))
               || (sd < 0.) ) 
              {
                sd = -b + d ; // 2nd root
              }
              if( sd > halfRmaxTolerance )  
              {
                if( sd < stheta )
                {
                  stheta    = sd;
                  sidetheta = kETheta;
                }
              }
            }
            else // sTheta+fDTheta > pi/2, concave surface, no normal
            {
              sd = -b - d;         // First root

              if ( ((std::fabs(sd) < halfRmaxTolerance) && (t2 >= 0.))
                || (sd < 0.)
                || ( (sd > 0.) && (p.z() + sd*v.z() > halfRmaxTolerance) ) )
              {
                sd = -b + d ; // 2nd root
              }
              if ( ( sd>halfRmaxTolerance )
                && ( p.z()+sd*v.z() <= halfRmaxTolerance ) )
              {
                if( sd < stheta )
                {
                  stheta    = sd;
                  sidetheta = kETheta;
                }
              }            
            }
          }
        }
      }
    }

  } // end theta intersections

  // Phi Intersection
    
  if ( !fFullPhiSphere )
  {
    if ( p.x() || p.y() ) // Check if on z axis (rho not needed later)
    {
      // pDist -ve when inside

      pDistS=p.x()*sinSPhi-p.y()*cosSPhi;
      pDistE=-p.x()*sinEPhi+p.y()*cosEPhi;

      // Comp -ve when in direction of outwards normal

      compS   = -sinSPhi*v.x()+cosSPhi*v.y() ;
      compE   =  sinEPhi*v.x()-cosEPhi*v.y() ;
      sidephi = kNull ;

      if ( (pDistS <= 0) && (pDistE <= 0) )
      {
        // Inside both phi *full* planes

        if ( compS < 0 )
        {
          sphi = pDistS/compS ;
          xi   = p.x()+sphi*v.x() ;
          yi   = p.y()+sphi*v.y() ;

          // Check intersection with correct half-plane (if not -> no intersect)
          //
          if( (std::fabs(xi)<=kCarTolerance) && (std::fabs(yi)<=kCarTolerance) )
          {
            vphi = std::atan2(v.y(),v.x());
            sidephi = kSPhi;
            if ( ( (fSPhi-halfAngTolerance) <= vphi)
              && ( (ePhi+halfAngTolerance)  >= vphi) )
            {
              sphi = kInfinity;
            }
          }
          else if ( ( yi*cosCPhi - xi*sinCPhi ) >= 0 )
          {
            sphi=kInfinity;
          }
          else
          {
            sidephi = kSPhi ;
            if ( pDistS > -halfCarTolerance)  { sphi = 0; } // Leave by sphi 
          }
        }
        else  { sphi = kInfinity; }

        if ( compE < 0 )
        {
          sphi2=pDistE/compE ;
          if (sphi2 < sphi) // Only check further if < starting phi intersection
          {
            xi = p.x()+sphi2*v.x() ;
            yi = p.y()+sphi2*v.y() ;

            // Check intersection with correct half-plane
            //
            if ( (std::fabs(xi)<=kCarTolerance)
              && (std::fabs(yi)<=kCarTolerance))
            {
              // Leaving via ending phi
              //
              vphi = std::atan2(v.y(),v.x()) ;
               
              if( !((fSPhi-halfAngTolerance <= vphi)
                  &&(fSPhi+fDPhi+halfAngTolerance >= vphi)) )
              { 
                sidephi = kEPhi;
                if ( pDistE <= -halfCarTolerance )  { sphi = sphi2; }
                else                                { sphi = 0.0;   }
              }
            }
            else if ((yi*cosCPhi-xi*sinCPhi)>=0) // Leaving via ending phi
            {
              sidephi = kEPhi ;
              if ( pDistE <= -halfCarTolerance )
              {
                sphi=sphi2;
              }
              else 
              {
                sphi = 0 ;
              }
            }
          }
        }        
      }
      else if ((pDistS >= 0) && (pDistE >= 0)) // Outside both *full* phi planes
      {
        if ( pDistS <= pDistE )
        {
          sidephi = kSPhi ;
        }
        else
        {
          sidephi = kEPhi ;
        }
        if ( fDPhi > pi )
        {
          if ( (compS < 0) && (compE < 0) )  { sphi = 0; }
          else                               { sphi = kInfinity; }
        }
        else
        {
          // if towards both >=0 then once inside (after error)
          // will remain inside

          if ( (compS >= 0) && (compE >= 0) ) { sphi = kInfinity; }
          else                                { sphi = 0; }
        }    
      }
      else if ( (pDistS > 0) && (pDistE < 0) )
      {
        // Outside full starting plane, inside full ending plane

        if ( fDPhi > pi )
        {
          if ( compE < 0 )
          {
            sphi = pDistE/compE ;
            xi   = p.x() + sphi*v.x() ;
            yi   = p.y() + sphi*v.y() ;

            // Check intersection in correct half-plane
            // (if not -> not leaving phi extent)
            //
            if( (std::fabs(xi)<=kCarTolerance)&&(std::fabs(yi)<=kCarTolerance) )
            {
              vphi = std::atan2(v.y(),v.x());
              sidephi = kSPhi;
              if ( ( (fSPhi-halfAngTolerance) <= vphi)
                && ( (ePhi+halfAngTolerance)  >= vphi) )
              {
                sphi = kInfinity;
              }
            }
            else if ( ( yi*cosCPhi - xi*sinCPhi ) <= 0 )
            {
              sphi = kInfinity ;
            }
            else // Leaving via Ending phi
            {
              sidephi = kEPhi ;
              if ( pDistE > -halfCarTolerance )  { sphi = 0.; }
            }
          }
          else
          {
            sphi = kInfinity ;
          }
        }
        else
        {
          if ( compS >= 0 )
          {
            if ( compE < 0 )
            {            
              sphi = pDistE/compE ;
              xi   = p.x() + sphi*v.x() ;
              yi   = p.y() + sphi*v.y() ;

              // Check intersection in correct half-plane
              // (if not -> remain in extent)
              //
              if( (std::fabs(xi)<=kCarTolerance)
               && (std::fabs(yi)<=kCarTolerance) )
              {
                vphi = std::atan2(v.y(),v.x());
                sidephi = kSPhi;
                if ( ( (fSPhi-halfAngTolerance) <= vphi)
                  && ( (ePhi+halfAngTolerance)  >= vphi) )
                {
                  sphi = kInfinity;
                }
              }
              else if ( ( yi*cosCPhi - xi*sinCPhi) <= 0 )
              {
                sphi=kInfinity;
              }
              else // otherwise leaving via Ending phi
              {
                sidephi = kEPhi ;
              }
            }
            else sphi=kInfinity;
          }
          else // leaving immediately by starting phi
          {
            sidephi = kSPhi ;
            sphi    = 0 ;
          }
        }
      }
      else
      {
        // Must be pDistS < 0 && pDistE > 0
        // Inside full starting plane, outside full ending plane

        if ( fDPhi > pi )
        {
          if ( compS < 0 )
          {
            sphi=pDistS/compS;
            xi=p.x()+sphi*v.x();
            yi=p.y()+sphi*v.y();
            
            // Check intersection in correct half-plane
            // (if not -> not leaving phi extent)
            //
            if( (std::fabs(xi)<=kCarTolerance)&&(std::fabs(yi)<=kCarTolerance) )
            {
              vphi = std::atan2(v.y(),v.x()) ;
              sidephi = kSPhi;
              if ( ( (fSPhi-halfAngTolerance) <= vphi)
                && ( (ePhi+halfAngTolerance)  >= vphi) )
              {
              sphi = kInfinity;
              }
            }
            else if ( ( yi*cosCPhi - xi*sinCPhi ) >= 0 )
            {
              sphi = kInfinity ;
            }
            else  // Leaving via Starting phi
            {
              sidephi = kSPhi ;
              if ( pDistS > -halfCarTolerance )  { sphi = 0; }
            }
          }
          else
          {
            sphi = kInfinity ;
          }
        }
        else
        {
          if ( compE >= 0 )
          {
            if ( compS < 0 )
            {
              sphi = pDistS/compS ;
              xi   = p.x()+sphi*v.x() ;
              yi   = p.y()+sphi*v.y() ;
              
              // Check intersection in correct half-plane
              // (if not -> remain in extent)
              //
              if( (std::fabs(xi)<=kCarTolerance)
               && (std::fabs(yi)<=kCarTolerance))
              {
                vphi = std::atan2(v.y(),v.x()) ;
                sidephi = kSPhi;
                if ( ( (fSPhi-halfAngTolerance) <= vphi)
                  && ( (ePhi+halfAngTolerance)  >= vphi) )
                {
                  sphi = kInfinity;
                }
              }
              else if ( ( yi*cosCPhi - xi*sinCPhi ) >= 0 )
              {
                sphi = kInfinity ;
              }
              else // otherwise leaving via Starting phi
              {
                sidephi = kSPhi ;
              }
            }
            else
            {
              sphi = kInfinity ;
            }
          }
          else // leaving immediately by ending
          {
            sidephi = kEPhi ;
            sphi    = 0     ;
          }
        }
      }      
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if ( v.x() || v.y() )
      {
        vphi = std::atan2(v.y(),v.x()) ;
        if ((fSPhi-halfAngTolerance < vphi) && (vphi < ePhi+halfAngTolerance))
        {
          sphi = kInfinity;
        }
        else
        {
          sidephi = kSPhi ; // arbitrary 
          sphi    = 0     ;
        }
      }
      else  // travel along z - no phi intersection
      {
        sphi = kInfinity ;
      }
    }
    if ( sphi < snxt )  // Order intersecttions
    {
      snxt = sphi ;
      side = sidephi ;
    }
  }
  if (stheta < snxt ) // Order intersections
  {
    snxt = stheta ;
    side = sidetheta ;
  }

  if (calcNorm)    // Output switch operator
  {
    switch( side )
    {
      case kRMax:
        xi=p.x()+snxt*v.x();
        yi=p.y()+snxt*v.y();
        zi=p.z()+snxt*v.z();
        *n=G4ThreeVector(xi/fRmax,yi/fRmax,zi/fRmax);
        *validNorm=true;
        break;

      case kRMin:
        *validNorm=false;  // Rmin is concave
        break;

      case kSPhi:
        if ( fDPhi <= pi )     // Normal to Phi-
        {
          *n=G4ThreeVector(sinSPhi,-cosSPhi,0);
          *validNorm=true;
        }
        else  { *validNorm=false; }
        break ;

      case kEPhi:
        if ( fDPhi <= pi )      // Normal to Phi+
        {
          *n=G4ThreeVector(-sinEPhi,cosEPhi,0);
          *validNorm=true;
        }
        else  { *validNorm=false; }
        break;

      case kSTheta:
        if( fSTheta == halfpi )
        {
          *n=G4ThreeVector(0.,0.,1.);
          *validNorm=true;
        }
        else if ( fSTheta > halfpi )
        {
          xi = p.x() + snxt*v.x();
          yi = p.y() + snxt*v.y();
          rho2=xi*xi+yi*yi;
          if (rho2)
          { 
            rhoSecTheta = std::sqrt(rho2*(1+tanSTheta2));
            *n = G4ThreeVector( xi/rhoSecTheta, yi/rhoSecTheta,
                               -tanSTheta/std::sqrt(1+tanSTheta2));
          }
          else
          {
            *n = G4ThreeVector(0.,0.,1.);
          }
          *validNorm=true;
        }
        else  { *validNorm=false; }  // Concave STheta cone
        break;

      case kETheta:
        if( eTheta == halfpi )
        {
          *n         = G4ThreeVector(0.,0.,-1.);
          *validNorm = true;
        }
        else if ( eTheta < halfpi )
        {
          xi=p.x()+snxt*v.x();
          yi=p.y()+snxt*v.y();
          rho2=xi*xi+yi*yi;
          if (rho2)
          { 
            rhoSecTheta = std::sqrt(rho2*(1+tanETheta2));
            *n = G4ThreeVector( xi/rhoSecTheta, yi/rhoSecTheta,
                               -tanETheta/std::sqrt(1+tanETheta2) );
          }
          else
          {
            *n = G4ThreeVector(0.,0.,-1.);
          }
          *validNorm=true;
        }
        else  { *validNorm=false; }   // Concave ETheta cone
        break;

      default:
        G4cout << G4endl;
        DumpInfo();
        std::ostringstream message;
        G4int oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x() = "   << p.x()/mm << " mm" << G4endl
                << "p.y() = "   << p.y()/mm << " mm" << G4endl
                << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
                << "Direction:" << G4endl << G4endl
                << "v.x() = "   << v.x() << G4endl
                << "v.y() = "   << v.y() << G4endl
                << "v.z() = "   << v.z() << G4endl << G4endl
                << "Proposed distance :" << G4endl << G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl;
        message.precision(oldprc);
        G4Exception("G4Sphere::DistanceToOut(p,v,..)",
                    "GeomSolids1002", JustWarning, message);
        break;
    }
  }
  if (snxt == kInfinity)
  {
    G4cout << G4endl;
    DumpInfo();
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Logic error: snxt = kInfinity  ???" << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Rp = "<< std::sqrt( p.x()*p.x()+p.y()*p.y()+p.z()*p.z() )/mm
            << " mm" << G4endl << G4endl
            << "Direction:" << G4endl << G4endl
            << "v.x() = "   << v.x() << G4endl
            << "v.y() = "   << v.y() << G4endl
            << "v.z() = "   << v.z() << G4endl << G4endl
            << "Proposed distance :" << G4endl << G4endl
            << "snxt = "    << snxt/mm << " mm" << G4endl;
    message.precision(oldprc);
    G4Exception("G4Sphere::DistanceToOut(p,v,..)",
                "GeomSolids1002", JustWarning, message);
  }

  return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Sphere::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0,safeRMin,safeRMax,safePhi,safeTheta;
  G4double rho2,rds,rho;
  G4double pTheta,dTheta1 = kInfinity,dTheta2 = kInfinity;
  rho2=p.x()*p.x()+p.y()*p.y();
  rds=std::sqrt(rho2+p.z()*p.z());
  rho=std::sqrt(rho2);

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
     G4int old_prc = G4cout.precision(16);
     G4cout << G4endl;
     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout.precision(old_prc) ;
     G4Exception("G4Sphere::DistanceToOut(p)",
                 "GeomSolids1002", JustWarning, "Point p is outside !?" );
  }
#endif

  // Distance to r shells
  //
  safeRMax = fRmax-rds;
  safe = safeRMax;  
  if (fRmin)
  {
     safeRMin = rds-fRmin;
     safe = std::min( safeRMin, safeRMax ); 
  }

  // Distance to phi extent
  //
  if ( !fFullPhiSphere )
  {
     if (rho>0.0)
     {
        if ((p.y()*cosCPhi-p.x()*sinCPhi)<=0)
        {
           safePhi=-(p.x()*sinSPhi-p.y()*cosSPhi);
        }
        else
        {
           safePhi=(p.x()*sinEPhi-p.y()*cosEPhi);
        }
     }
     else
     {
        safePhi = 0.0;  // Distance to both Phi surfaces (extended)
     }
     // Both cases above can be improved - in case fRMin > 0.0
     //  although it may be costlier (good for precise, not fast version)
     
     safe= std::min(safe, safePhi);
  }

  // Distance to Theta extent
  //
  if ( !fFullThetaSphere )
  {
    if( rds > 0.0 )
    {
       pTheta=std::acos(p.z()/rds);
       if (pTheta<0) { pTheta+=pi; }
       if(fSTheta>0.)
       { dTheta1=pTheta-fSTheta;}
       if(eTheta<pi)
       { dTheta2=eTheta-pTheta;}
      
       safeTheta=rds*std::sin(std::min(dTheta1, dTheta2) );
    }
    else
    {
       safeTheta= 0.0;
         // An improvement will be to return negative answer if outside (TODO)
    }
    safe = std::min( safe, safeTheta );
  }

  if (safe<0.0) { safe=0; }
    // An improvement to return negative answer if outside (TODO)
  
  return safe;
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType

G4GeometryType G4Sphere::GetEntityType() const
{
  return G4String("G4Sphere");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Sphere::Clone() const
{
  return new G4Sphere(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Sphere::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Sphere\n"
     << " Parameters: \n"
     << "    inner radius: " << fRmin/mm << " mm \n"
     << "    outer radius: " << fRmax/mm << " mm \n"
     << "    starting phi of segment  : " << fSPhi/degree << " degrees \n"
     << "    delta phi of segment     : " << fDPhi/degree << " degrees \n"
     << "    starting theta of segment: " << fSTheta/degree << " degrees \n"
     << "    delta theta of segment   : " << fDTheta/degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

////////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Sphere::GetPointOnSurface() const
{
  G4double zRand, aOne, aTwo, aThr, aFou, aFiv, chose, phi, sinphi, cosphi;
  G4double height1, height2, slant1, slant2, costheta, sintheta, rRand;

  height1 = (fRmax-fRmin)*cosSTheta;
  height2 = (fRmax-fRmin)*cosETheta;
  slant1  = std::sqrt(sqr((fRmax - fRmin)*sinSTheta) + height1*height1);
  slant2  = std::sqrt(sqr((fRmax - fRmin)*sinETheta) + height2*height2);
  rRand   = GetRadiusInRing(fRmin,fRmax);
  
  aOne = fRmax*fRmax*fDPhi*(cosSTheta-cosETheta);
  aTwo = fRmin*fRmin*fDPhi*(cosSTheta-cosETheta);
  aThr = fDPhi*((fRmax + fRmin)*sinSTheta)*slant1;
  aFou = fDPhi*((fRmax + fRmin)*sinETheta)*slant2;
  aFiv = 0.5*fDTheta*(fRmax*fRmax-fRmin*fRmin);
  
  phi = G4RandFlat::shoot(fSPhi, ePhi); 
  cosphi = std::cos(phi); 
  sinphi = std::sin(phi);
  costheta = G4RandFlat::shoot(cosETheta,cosSTheta);
  sintheta = std::sqrt(1.-sqr(costheta));

  if(fFullPhiSphere) { aFiv = 0; }
  if(fSTheta == 0)   { aThr=0; }
  if(eTheta == pi) { aFou = 0; }
  if(fSTheta == halfpi) { aThr = pi*(fRmax*fRmax-fRmin*fRmin); }
  if(eTheta == halfpi)  { aFou = pi*(fRmax*fRmax-fRmin*fRmin); }

  chose = G4RandFlat::shoot(0.,aOne+aTwo+aThr+aFou+2.*aFiv);
  if( (chose>=0.) && (chose<aOne) )
  {
    return G4ThreeVector(fRmax*sintheta*cosphi,
                         fRmax*sintheta*sinphi, fRmax*costheta);
  }
  else if( (chose>=aOne) && (chose<aOne+aTwo) )
  {
    return G4ThreeVector(fRmin*sintheta*cosphi,
                         fRmin*sintheta*sinphi, fRmin*costheta);
  }
  else if( (chose>=aOne+aTwo) && (chose<aOne+aTwo+aThr) )
  {
    if (fSTheta != halfpi)
    {
      zRand = G4RandFlat::shoot(fRmin*cosSTheta,fRmax*cosSTheta);
      return G4ThreeVector(tanSTheta*zRand*cosphi,
                           tanSTheta*zRand*sinphi,zRand);
    }
    else
    {
      return G4ThreeVector(rRand*cosphi, rRand*sinphi, 0.);
    }    
  }
  else if( (chose>=aOne+aTwo+aThr) && (chose<aOne+aTwo+aThr+aFou) )
  {
    if(eTheta != halfpi)
    {
      zRand = G4RandFlat::shoot(fRmin*cosETheta, fRmax*cosETheta);
      return G4ThreeVector  (tanETheta*zRand*cosphi,
                             tanETheta*zRand*sinphi,zRand);
    }
    else
    {
      return G4ThreeVector(rRand*cosphi, rRand*sinphi, 0.);
    }
  }
  else if( (chose>=aOne+aTwo+aThr+aFou) && (chose<aOne+aTwo+aThr+aFou+aFiv) )
  {
    return G4ThreeVector(rRand*sintheta*cosSPhi,
                         rRand*sintheta*sinSPhi,rRand*costheta);
  }
  else
  {
    return G4ThreeVector(rRand*sintheta*cosEPhi,
                         rRand*sintheta*sinEPhi,rRand*costheta);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// GetSurfaceArea

G4double G4Sphere::GetSurfaceArea()
{
  if(fSurfaceArea != 0.) {;}
  else
  {   
    G4double Rsq=fRmax*fRmax;
    G4double rsq=fRmin*fRmin;
         
    fSurfaceArea = fDPhi*(rsq+Rsq)*(cosSTheta - cosETheta);
    if(!fFullPhiSphere)
    {
      fSurfaceArea = fSurfaceArea + fDTheta*(Rsq-rsq);
    }
    if(fSTheta >0)
    {
      G4double acos1=std::acos( std::pow(sinSTheta,2) * std::cos(fDPhi)
                              + std::pow(cosSTheta,2));
      if(fDPhi>pi)
      { 
        fSurfaceArea = fSurfaceArea + 0.5*(Rsq-rsq)*(twopi-acos1);
      }
      else
      {
        fSurfaceArea = fSurfaceArea + 0.5*(Rsq-rsq)*acos1;
      }
    }
    if(eTheta < pi)
    {
      G4double acos2=std::acos( std::pow(sinETheta,2) * std::cos(fDPhi)
                              + std::pow(cosETheta,2));
      if(fDPhi>pi)
      { 
        fSurfaceArea = fSurfaceArea + 0.5*(Rsq-rsq)*(twopi-acos2);
      }
      else
      {
        fSurfaceArea = fSurfaceArea + 0.5*(Rsq-rsq)*acos2;
      }
    }
  }
  return fSurfaceArea;
}

/////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

G4VisExtent G4Sphere::GetExtent() const
{
  return G4VisExtent(-fRmax, fRmax,-fRmax, fRmax,-fRmax, fRmax );
}


void G4Sphere::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Sphere::CreatePolyhedron () const
{
  return new G4PolyhedronSphere (fRmin, fRmax, fSPhi, fDPhi, fSTheta, fDTheta);
}

#endif
