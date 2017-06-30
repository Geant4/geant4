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
// $Id: G4CutTubs.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// class G4CutTubs
//
// History:
//
// 30.10.16 E.Tcherniaev - reimplemented CalculateExtent(),
//                       removed CreateRotatedVetices()
// 05.04.12 M.Kelsey   - GetPointOnSurface() throw flat in sqrt(r)
// 01.06.11 T.Nikitina - Derived from G4Tubs
//
/////////////////////////////////////////////////////////////////////////

#include "G4CutTubs.hh"

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4CutTubs::G4CutTubs( const G4String &pName,
                      G4double pRMin, G4double pRMax,
                      G4double pDz,
                      G4double pSPhi, G4double pDPhi,
                      G4ThreeVector pLowNorm,G4ThreeVector pHighNorm )
  : G4OTubs(pName, pRMin, pRMax, pDz, pSPhi, pDPhi),
    fPhiFullCutTube(true)
{
  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  halfCarTolerance = kCarTolerance*0.5;
  halfRadTolerance = kRadTolerance*0.5;
  halfAngTolerance = kAngTolerance*0.5;

  // Check on Cutted Planes Normals
  // If there is NO CUT, propose to use G4Tubs instead
  //
  if(pDPhi<twopi)  { fPhiFullCutTube=false; }
 
  if ( ( !pLowNorm.x()) && ( !pLowNorm.y())
    && ( !pHighNorm.x()) && (!pHighNorm.y()) )
  {
    std::ostringstream message;
    message << "Inexisting Low/High Normal to Z plane or Parallel to Z."
            << G4endl
            << "Normals to Z plane are (" << pLowNorm <<" and "
            << pHighNorm << ") in solid: " << GetName();
    G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids1001",
                JustWarning, message, "Should use G4Tubs!");
  }

  // If Normal is (0,0,0),means parallel to R, give it value of (0,0,+/-1)
  // 
  if (pLowNorm.mag2() == 0.)  { pLowNorm.setZ(-1.); }
  if (pHighNorm.mag2()== 0.)  { pHighNorm.setZ(1.); }

  // Given Normals to Cut Planes have to be an unit vectors. 
  // Normalize if it is needed.
  //
  if (pLowNorm.mag2() != 1.)  { pLowNorm  = pLowNorm.unit();  }
  if (pHighNorm.mag2()!= 1.)  { pHighNorm = pHighNorm.unit(); }

  // Normals to cutted planes have to point outside Solid
  //
  if( (pLowNorm.mag2() != 0.) && (pHighNorm.mag2()!= 0. ) )
  {
    if( ( pLowNorm.z()>= 0. ) || ( pHighNorm.z() <= 0.))
    {
      std::ostringstream message;
      message << "Invalid Low or High Normal to Z plane; "
                 "has to point outside Solid." << G4endl
              << "Invalid Norm to Z plane (" << pLowNorm << " or  "
              << pHighNorm << ") in solid: " << GetName();
      G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids0002",
                  FatalException, message);
    }
  }
  fLowNorm  = pLowNorm;
  fHighNorm = pHighNorm;

  // Check Intersection of cut planes. They MUST NOT Intersect
  //
  // This check has been disabled as too strict.
  // See problem report #1887
  //
  // if(IsCrossingCutPlanes())
  // {
  //   std::ostringstream message;
  //   message << "Invalid Low or High Normal to Z plane; "
  //           << "Crossing Cutted Planes." << G4endl
  //           << "Invalid Norm to Z plane (" << pLowNorm << " and "
  //           << pHighNorm << ") in solid: " << GetName();
  //   G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids0002",
  //               FatalException, message);
  // }
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4CutTubs::G4CutTubs( __void__& a )
  : G4OTubs(a), fLowNorm(G4ThreeVector()),
    fHighNorm(G4ThreeVector()), fPhiFullCutTube(false),
    halfCarTolerance(0.), halfRadTolerance(0.), halfAngTolerance(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4CutTubs::~G4CutTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4CutTubs::G4CutTubs(const G4CutTubs& rhs)
  : G4OTubs(rhs), fLowNorm(rhs.fLowNorm), fHighNorm(rhs.fHighNorm),
    fPhiFullCutTube(rhs.fPhiFullCutTube),
    halfCarTolerance(rhs.halfCarTolerance),
    halfRadTolerance(rhs.halfRadTolerance),
    halfAngTolerance(rhs.halfAngTolerance)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4CutTubs& G4CutTubs::operator = (const G4CutTubs& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4OTubs::operator=(rhs);

   // Copy data
   //
   fLowNorm = rhs.fLowNorm; fHighNorm = rhs.fHighNorm;
   fPhiFullCutTube = rhs.fPhiFullCutTube;
   halfCarTolerance = rhs.halfCarTolerance;
   halfRadTolerance = rhs.halfRadTolerance;
   halfAngTolerance = rhs.halfAngTolerance;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4CutTubs::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double rmin = GetInnerRadius();
  G4double rmax = GetOuterRadius();
  G4double dz   = GetZHalfLength();

  G4ThreeVector norm;
  G4double xynorm, znorm;

  // get zmin
  norm = GetLowNorm();
  xynorm = std::sqrt(norm.x()*norm.x()+norm.y()*norm.y());
  znorm  = std::abs(norm.z());
  G4double zmin = -(dz + rmax*xynorm/znorm); 

  // get zmax
  norm = GetHighNorm();
  xynorm = std::sqrt(norm.x()*norm.x()+norm.y()*norm.y());
  znorm  = std::abs(norm.z());
  G4double zmax = dz + rmax*xynorm/znorm; 

  // Find bounding box
  //
  if (GetDeltaPhiAngle() < twopi)
  {
    G4TwoVector vmin,vmax;
    G4GeomTools::DiskExtent(rmin,rmax,
                            GetSinStartPhi(),GetCosStartPhi(),
                            GetSinEndPhi(),GetCosEndPhi(),
                            vmin,vmax);
    pMin.set(vmin.x(),vmin.y(), zmin);
    pMax.set(vmax.x(),vmax.y(), zmax);
  }
  else
  {
    pMin.set(-rmax,-rmax, zmin);
    pMax.set( rmax, rmax, zmax);
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
    G4Exception("G4CutTubs::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4CutTubs::CalculateExtent( const EAxis              pAxis,
                                   const G4VoxelLimits&     pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double&          pMin, 
                                         G4double&          pMax    ) const
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
  G4double dphi = GetDeltaPhiAngle();
  G4double zmin = bmin.z();
  G4double zmax = bmax.z();

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
      baseA[k].set(rext*cosCur,rext*sinCur,zmin);
      baseB[k].set(rext*cosCur,rext*sinCur,zmax);

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
    pols[0][0].set(rmin*cosStart,rmin*sinStart,zmax);
    pols[0][1].set(rmin*cosStart,rmin*sinStart,zmin);
    pols[0][2].set(rmax*cosStart,rmax*sinStart,zmin);
    pols[0][3].set(rmax*cosStart,rmax*sinStart,zmax);
    for (G4int k=1; k<ksteps+1; ++k)
    {
      pols[k][0].set(rmin*cosCur,rmin*sinCur,zmax);
      pols[k][1].set(rmin*cosCur,rmin*sinCur,zmin);
      pols[k][2].set(rext*cosCur,rext*sinCur,zmin);
      pols[k][3].set(rext*cosCur,rext*sinCur,zmax);

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    pols[ksteps+1][0].set(rmin*cosEnd,rmin*sinEnd,zmax);
    pols[ksteps+1][1].set(rmin*cosEnd,rmin*sinEnd,zmin);
    pols[ksteps+1][2].set(rmax*cosEnd,rmax*sinEnd,zmin);
    pols[ksteps+1][3].set(rmax*cosEnd,rmax*sinEnd,zmax);

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
// Return whether point inside/outside/on surface

EInside G4CutTubs::Inside( const G4ThreeVector& p ) const
{
  G4ThreeVector vZ = G4ThreeVector(0,0,fDz);
  EInside in = kInside;

  // Check the lower cut plane
  //
  G4double zinLow =(p+vZ).dot(fLowNorm);
  if (zinLow > halfCarTolerance)  { return kOutside; }

  // Check the higher cut plane
  //
  G4double zinHigh = (p-vZ).dot(fHighNorm);
  if (zinHigh > halfCarTolerance)  { return kOutside; }

  // Check radius
  //
  G4double r2 = p.x()*p.x() + p.y()*p.y() ;

  G4double tolRMin = fRMin - halfRadTolerance;
  G4double tolRMax = fRMax + halfRadTolerance;
  if ( tolRMin < 0 )  { tolRMin = 0; }

  if (r2 < tolRMin*tolRMin || r2 > tolRMax*tolRMax) { return kOutside; }

  // Check Phi cut
  //
  if(!fPhiFullCutTube)
  {
    if ((tolRMin == 0) && (std::fabs(p.x()) <= halfCarTolerance)
                       && (std::fabs(p.y()) <= halfCarTolerance))
    {
      return kSurface;
    }

    G4double phi0 = std::atan2(p.y(),p.x());
    G4double phi1 = phi0 - twopi;
    G4double phi2 = phi0 + twopi;

    in = kOutside;
    G4double sphi = fSPhi - halfAngTolerance;
    G4double ephi = sphi + fDPhi + kAngTolerance;
    if ((phi0  >= sphi && phi0  <= ephi) ||
        (phi1  >= sphi && phi1  <= ephi) ||
        (phi2  >= sphi && phi2  <= ephi)) in = kSurface;
    if (in == kOutside)  { return kOutside; }

    sphi += kAngTolerance;
    ephi -= kAngTolerance;
    if ((phi0  >= sphi && phi0  <= ephi) ||
        (phi1  >= sphi && phi1  <= ephi) ||
        (phi2  >= sphi && phi2  <= ephi)) in = kInside;
    if (in == kSurface)  { return kSurface; }
  }

  // Check on the Surface for Z
  //
  if ((zinLow >= -halfCarTolerance) || (zinHigh >= -halfCarTolerance))
  {
    return kSurface;
  }

  // Check on the Surface for R
  //
  if (fRMin) { tolRMin = fRMin + halfRadTolerance; }
  else       { tolRMin = 0; }
  tolRMax = fRMax - halfRadTolerance;
  if (((r2 <= tolRMin*tolRMin) || (r2 >= tolRMax*tolRMax)) &&
       (r2 >= halfRadTolerance*halfRadTolerance))
  {
    return kSurface;
  }

  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4CutTubs::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4int noSurfaces = 0;
  G4double rho, pPhi;
  G4double distZLow,distZHigh, distRMin, distRMax;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  G4ThreeVector norm, sumnorm(0.,0.,0.);
  G4ThreeVector nZ = G4ThreeVector(0, 0, 1.0);
  G4ThreeVector nR, nPs, nPe;

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y());

  distRMin = std::fabs(rho - fRMin);
  distRMax = std::fabs(rho - fRMax);

  // dist to Low Cut
  //
  distZLow =std::fabs((p+vZ).dot(fLowNorm));
 
  // dist to High Cut
  //
  distZHigh = std::fabs((p-vZ).dot(fHighNorm));

  if (!fPhiFullCutTube)    // Protected against (0,0,z) 
  {
    if ( rho > halfCarTolerance )
    {
      pPhi = std::atan2(p.y(),p.x());
    
      if(pPhi  < fSPhi- halfCarTolerance)           { pPhi += twopi; }
      else if(pPhi > fSPhi+fDPhi+ halfCarTolerance) { pPhi -= twopi; }

      distSPhi = std::fabs(pPhi - fSPhi);       
      distEPhi = std::fabs(pPhi - fSPhi - fDPhi); 
    }
    else if( !fRMin )
    {
      distSPhi = 0.; 
      distEPhi = 0.; 
    }
    nPs = G4ThreeVector(std::sin(fSPhi),-std::cos(fSPhi),0);
    nPe = G4ThreeVector(-std::sin(fSPhi+fDPhi),std::cos(fSPhi+fDPhi),0);
  }
  if ( rho > halfCarTolerance ) { nR = G4ThreeVector(p.x()/rho,p.y()/rho,0); }

  if( distRMax <= halfCarTolerance ) 
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  if( fRMin && (distRMin <= halfCarTolerance) )
  {
    noSurfaces ++;
    sumnorm -= nR;
  }
  if( fDPhi < twopi )   
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
  if (distZLow <= halfCarTolerance)  
  {
    noSurfaces ++;
    sumnorm += fLowNorm;
  }
  if (distZHigh <= halfCarTolerance)  
  {
    noSurfaces ++;
    sumnorm += fHighNorm;
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4CutTubs::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" );
    G4int oldprc = G4cout.precision(20);
    G4cout<< "G4CutTubs::SN ( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "
          << G4endl << G4endl;
    G4cout.precision(oldprc) ;
#endif 
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = sumnorm.unit(); }

  return norm;
}

/////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4CutTubs::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  ENorm side ;
  G4ThreeVector norm ;
  G4double rho, phi ;
  G4double distZLow,distZHigh,distZ;
  G4double distRMin, distRMax, distSPhi, distEPhi, distMin ;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

  distRMin = std::fabs(rho - fRMin) ;
  distRMax = std::fabs(rho - fRMax) ;

  //dist to Low Cut
  //
  distZLow =std::fabs((p+vZ).dot(fLowNorm));

  //dist to High Cut
  //
  distZHigh = std::fabs((p-vZ).dot(fHighNorm));
  distZ=std::min(distZLow,distZHigh);

  if (distRMin < distRMax) // First minimum
  {
    if ( distZ < distRMin )
    {
       distMin = distZ ;
       side    = kNZ ;
    }
    else
    {
      distMin = distRMin ;
      side    = kNRMin   ;
    }
  }
  else
  {
    if ( distZ < distRMax )
    {
      distMin = distZ ;
      side    = kNZ   ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax   ;
    }
  }   
  if (!fPhiFullCutTube  &&  rho ) // Protected against (0,0,z) 
  {
    phi = std::atan2(p.y(),p.x()) ;

    if ( phi < 0 )  { phi += twopi; }

    if ( fSPhi < 0 )
    {
      distSPhi = std::fabs(phi - (fSPhi + twopi))*rho ;
    }
    else
    {
      distSPhi = std::fabs(phi - fSPhi)*rho ;
    }
    distEPhi = std::fabs(phi - fSPhi - fDPhi)*rho ;
                                      
    if (distSPhi < distEPhi) // Find new minimum
    {
      if ( distSPhi < distMin )
      {
        side = kNSPhi ;
      }
    }
    else
    {
      if ( distEPhi < distMin )
      {
        side = kNEPhi ;
      }
    }
  }    
  switch ( side )
  {
    case kNRMin : // Inner radius
    {                      
      norm = G4ThreeVector(-p.x()/rho, -p.y()/rho, 0) ;
      break ;
    }
    case kNRMax : // Outer radius
    {                  
      norm = G4ThreeVector(p.x()/rho, p.y()/rho, 0) ;
      break ;
    }
    case kNZ :    // + or - dz
    {                              
      if ( distZHigh > distZLow )  { norm = fHighNorm ; }
      else                         { norm = fLowNorm; }
      break ;
    }
    case kNSPhi:
    {
      norm = G4ThreeVector(std::sin(fSPhi), -std::cos(fSPhi), 0) ;
      break ;
    }
    case kNEPhi:
    {
      norm = G4ThreeVector(-std::sin(fSPhi+fDPhi), std::cos(fSPhi+fDPhi), 0) ;
      break;
    }
    default:      // Should never reach this case ...
    {
      DumpInfo();
      G4Exception("G4CutTubs::ApproxSurfaceNormal()",
                  "GeomSolids1002", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break ;
    }    
  }                
  return norm;
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - 'if valid' implies tolerant checking of intersection points

G4double G4CutTubs::DistanceToIn( const G4ThreeVector& p,
                                  const G4ThreeVector& v  ) const
{
  G4double snxt = kInfinity ;      // snxt = default return value
  G4double tolORMin2, tolIRMax2 ;  // 'generous' radii squared
  G4double tolORMax2, tolIRMin2;
  const G4double dRmax = 100.*fRMax;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);
  
  // Intersection point variables
  //
  G4double Dist, sd=0, xi, yi, zi, rho2, inum, iden, cosPsi, Comp,calf ;
  G4double t1, t2, t3, b, c, d ;     // Quadratic solver variables 
  G4double distZLow,distZHigh;
  // Calculate tolerant rmin and rmax

  if (fRMin > kRadTolerance)
  {
    tolORMin2 = (fRMin - halfRadTolerance)*(fRMin - halfRadTolerance) ;
    tolIRMin2 = (fRMin + halfRadTolerance)*(fRMin + halfRadTolerance) ;
  }
  else
  {
    tolORMin2 = 0.0 ;
    tolIRMin2 = 0.0 ;
  }
  tolORMax2 = (fRMax + halfRadTolerance)*(fRMax + halfRadTolerance) ;
  tolIRMax2 = (fRMax - halfRadTolerance)*(fRMax - halfRadTolerance) ;

  // Intersection with ZCut surfaces

  // dist to Low Cut
  //
  distZLow =(p+vZ).dot(fLowNorm);

  // dist to High Cut
  //
  distZHigh = (p-vZ).dot(fHighNorm);

  if ( distZLow >= -halfCarTolerance )
  {
    calf = v.dot(fLowNorm);
    if (calf<0)
    {
      sd = -distZLow/calf;
      if(sd < 0.0)  { sd = 0.0; }

      xi   = p.x() + sd*v.x() ;                // Intersection coords
      yi   = p.y() + sd*v.y() ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
        if (!fPhiFullCutTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi*cosCPhi + yi*sinCPhi ;
          iden   = std::sqrt(rho2) ;
          cosPsi = inum/iden ;
          if (cosPsi >= cosHDPhiIT)  { return sd ; }
        }
        else
        {
          return sd ;
        }
      }
    }
    else
    {
      if ( sd<halfCarTolerance )
      {
        if(calf>=0) { sd=kInfinity; }
        return sd ;  // On/outside extent, and heading away
      }              // -> cannot intersect
    }
  }

  if(distZHigh >= -halfCarTolerance )
  {
    calf = v.dot(fHighNorm);
    if (calf<0)
    {
      sd = -distZHigh/calf;

      if(sd < 0.0)  { sd = 0.0; }

      xi   = p.x() + sd*v.x() ;                // Intersection coords
      yi   = p.y() + sd*v.y() ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
        if (!fPhiFullCutTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi*cosCPhi + yi*sinCPhi ;
          iden   = std::sqrt(rho2) ;
          cosPsi = inum/iden ;
          if (cosPsi >= cosHDPhiIT)  { return sd ; }
        }
        else
        {
          return sd ;
        }
      }
    }
    else
    {
      if ( sd<halfCarTolerance )
      { 
        if(calf>=0) { sd=kInfinity; }
        return sd ;  // On/outside extent, and heading away
      }              // -> cannot intersect
    }
  }

  // -> Can not intersect z surfaces
  //
  // Intersection with rmax (possible return) and rmin (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //            t1                t2                t3

  t1 = 1.0 - v.z()*v.z() ;
  t2 = p.x()*v.x() + p.y()*v.y() ;
  t3 = p.x()*p.x() + p.y()*p.y() ;
  if ( t1 > 0 )        // Check not || to z axis
  {
    b = t2/t1 ;
    c = t3 - fRMax*fRMax ;
    
    if ((t3 >= tolORMax2) && (t2<0))   // This also handles the tangent case
    {
      // Try outer cylinder intersection, c=(t3-fRMax*fRMax)/t1;

      c /= t1 ;
      d = b*b - c ;

      if (d >= 0)  // If real root
      {
        sd = c/(-b+std::sqrt(d));
        if (sd >= 0)  // If 'forwards'
        {
          if ( sd>dRmax ) // Avoid rounding errors due to precision issues on
          {               // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd-std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
          // Check z intersection
          //
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              // Z ok. Check phi intersection if reqd
              //
              if (fPhiFullCutTube)
              {
                return sd ;
              }
              else
              {
                xi     = p.x() + sd*v.x() ;
                yi     = p.y() + sd*v.y() ;
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/fRMax ;
                if (cosPsi >= cosHDPhiIT)  { return sd ; }
              }
            }  //  end if std::fabs(zi)
          }
        }    //  end if (sd>=0)
      }      //  end if (d>=0)
    }        //  end if (r>=fRMax)
    else 
    {
      // Inside outer radius :
      // check not inside, and heading through tubs (-> 0 to in)
      if ((t3 > tolIRMin2) && (t2 < 0)
       && (std::fabs(p.z()) <= std::fabs(GetCutZ(p))-halfCarTolerance ))
      {
        // Inside both radii, delta r -ve, inside z extent

        if (!fPhiFullCutTube)
        {
          inum   = p.x()*cosCPhi + p.y()*sinCPhi ;
          iden   = std::sqrt(t3) ;
          cosPsi = inum/iden ;
          if (cosPsi >= cosHDPhiIT)
          {
            // In the old version, the small negative tangent for the point
            // on surface was not taken in account, and returning 0.0 ...
            // New version: check the tangent for the point on surface and 
            // if no intersection, return kInfinity, if intersection instead
            // return sd.
            //
            c = t3-fRMax*fRMax; 
            if ( c<=0.0 )
            {
              return 0.0;
            }
            else
            {
              c = c/t1 ;
              d = b*b-c;
              if ( d>=0.0 )
              {
                snxt = c/(-b+std::sqrt(d)); // using safe solution
                                            // for quadratic equation 
                if ( snxt < halfCarTolerance ) { snxt=0; }
                return snxt ;
              }      
              else
              {
                return kInfinity;
              }
            }
          } 
        }
        else
        {   
          // In the old version, the small negative tangent for the point
          // on surface was not taken in account, and returning 0.0 ...
          // New version: check the tangent for the point on surface and 
          // if no intersection, return kInfinity, if intersection instead
          // return sd.
          //
          c = t3 - fRMax*fRMax; 
          if ( c<=0.0 )
          {
            return 0.0;
          }
          else
          {
            c = c/t1 ;
            d = b*b-c;
            if ( d>=0.0 )
            {
              snxt= c/(-b+std::sqrt(d)); // using safe solution
                                         // for quadratic equation 
              if ( snxt < halfCarTolerance ) { snxt=0; }
              return snxt ;
            }      
            else
            {
              return kInfinity;
            }
          }
        } // end if   (!fPhiFullCutTube)
      }   // end if   (t3>tolIRMin2)
    }     // end if   (Inside Outer Radius) 
      
    if ( fRMin )    // Try inner cylinder intersection
    {
      c = (t3 - fRMin*fRMin)/t1 ;
      d = b*b - c ;
      if ( d >= 0.0 )  // If real root
      {
        // Always want 2nd root - we are outside and know rmax Hit was bad
        // - If on surface of rmin also need farthest root
        
        sd =( b > 0. )? c/(-b - std::sqrt(d)) : (-b + std::sqrt(d));
        if (sd >= -10*halfCarTolerance)  // check forwards
        {
          // Check z intersection
          //
          if (sd < 0.0)  { sd = 0.0; }
          if (sd>dRmax) // Avoid rounding errors due to precision issues seen
          {             // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd-std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              // Z ok. Check phi
              //
              if ( fPhiFullCutTube )
              {
                return sd ; 
              }
              else
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/fRMin ;
                if (cosPsi >= cosHDPhiIT)
                {
                  // Good inner radius isect
                  // - but earlier phi isect still possible
                  //
                  snxt = sd ;
                }
              }
            }      //    end if std::fabs(zi)
          }
        }          //    end if (sd>=0)
      }            //    end if (d>=0)
    }              //    end if (fRMin)
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to kCarTolerance*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?
  //
  if ( !fPhiFullCutTube )
  {
    // First phi surface (Starting phi)
    //
    Comp = v.x()*sinSPhi - v.y()*cosSPhi ;
                
    if ( Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = (p.y()*cosSPhi - p.x()*sinSPhi) ;

      if ( Dist < halfCarTolerance )
      {
        sd = Dist/Comp ;

        if (sd < snxt)
        {
          if ( sd < 0 )  { sd = 0.0; }
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance) 
            {
              rho2 = xi*xi + yi*yi ;
              if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
                || ( (rho2 >  tolORMin2) && (rho2 <  tolIRMin2)
                  && ( v.y()*cosSPhi - v.x()*sinSPhi >  0 )
                  && ( v.x()*cosSPhi + v.y()*sinSPhi >= 0 )     )
                || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                  && (v.y()*cosSPhi - v.x()*sinSPhi > 0)
                  && (v.x()*cosSPhi + v.y()*sinSPhi < 0) )    )
              {
                // z and r intersections good
                // - check intersecting with correct half-plane
                //
                if ((yi*cosCPhi-xi*sinCPhi) <= halfCarTolerance) { snxt = sd; }
              }
            }   //two Z conditions
          }
        }
      }    
    }
      
    // Second phi surface (Ending phi)
    //
    Comp = -(v.x()*sinEPhi - v.y()*cosEPhi) ;
        
    if (Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = -(p.y()*cosEPhi - p.x()*sinEPhi) ;

      if ( Dist < halfCarTolerance )
      {
        sd = Dist/Comp ;

        if (sd < snxt)
        {
          if ( sd < 0 )  { sd = 0; }
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              xi   = p.x() + sd*v.x() ;
              yi   = p.y() + sd*v.y() ;
              rho2 = xi*xi + yi*yi ;
              if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
                  || ( (rho2 > tolORMin2)  && (rho2 < tolIRMin2)
                    && (v.x()*sinEPhi - v.y()*cosEPhi >  0)
                    && (v.x()*cosEPhi + v.y()*sinEPhi >= 0) )
                  || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                    && (v.x()*sinEPhi - v.y()*cosEPhi > 0)
                    && (v.x()*cosEPhi + v.y()*sinEPhi < 0) ) )
              {
                // z and r intersections good
                // - check intersecting with correct half-plane
                //
                if ( (yi*cosCPhi-xi*sinCPhi) >= -halfCarTolerance )
                {
                  snxt = sd;
                }
              }    //?? >=-halfCarTolerance
            }
          }  // two Z conditions
        }
      }
    }         //  Comp < 0
  }           //  !fPhiFullTube 
  if ( snxt<halfCarTolerance )  { snxt=0; }

  return snxt ;
}
 
//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points
//   Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4CutTubs::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safRMin,safRMax,safZLow,safZHigh,safePhi,safe,rho,cosPsi;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  // Distance to R
  //
  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

  safRMin = fRMin- rho ;
  safRMax = rho - fRMax ;

  // Distances to ZCut(Low/High)

  // Dist to Low Cut
  //
  safZLow = (p+vZ).dot(fLowNorm);

  // Dist to High Cut
  //
  safZHigh = (p-vZ).dot(fHighNorm);

  safe = std::max(safZLow,safZHigh);

  if ( safRMin > safe ) { safe = safRMin; }
  if ( safRMax> safe )  { safe = safRMax; }

  // Distance to Phi
  //
  if ( (!fPhiFullCutTube) && (rho) )
   {
     // Psi=angle from central phi to point
     //
     cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/rho ;
     
     if ( cosPsi < std::cos(fDPhi*0.5) )
     {
       // Point lies outside phi range
 
       if ( (p.y()*cosCPhi - p.x()*sinCPhi) <= 0 )
       {
         safePhi = std::fabs(p.x()*sinSPhi - p.y()*cosSPhi) ;
       }
       else
       {
         safePhi = std::fabs(p.x()*sinEPhi - p.y()*cosEPhi) ;
       }
       if ( safePhi > safe )  { safe = safePhi; }
     }
   }
   if ( safe < 0 )  { safe = 0; }

   return safe ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4CutTubs::DistanceToOut( const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm,
                                         G4bool *validNorm,
                                         G4ThreeVector *n    ) const
{  
  ESide side=kNull , sider=kNull, sidephi=kNull ;
  G4double snxt=kInfinity, srd=kInfinity,sz=kInfinity, sphi=kInfinity ;
  G4double deltaR, t1, t2, t3, b, c, d2, roMin2 ;
  G4double distZLow,distZHigh,calfH,calfL;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);
 
  // Vars for phi intersection:
  //
  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, vphi, roi2 ;
 
  // Z plane intersection
  // Distances to ZCut(Low/High)

  // dist to Low Cut
  //
  distZLow =(p+vZ).dot(fLowNorm);

  // dist to High Cut
  //
  distZHigh = (p-vZ).dot(fHighNorm);

  calfH = v.dot(fHighNorm);
  calfL = v.dot(fLowNorm);

  if (calfH > 0 )
  {
    if ( distZHigh < halfCarTolerance )
    {
      snxt = -distZHigh/calfH ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector(0,0,1) ;
        *validNorm = true ;
      }
      return snxt = 0 ;
    }
 }
  if ( calfL>0)
  {
   
    if ( distZLow < halfCarTolerance )
    {
      sz = -distZLow/calfL ;
      if(sz<snxt){
      snxt=sz;
      side = kMZ ;
      }
      
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector(0,0,-1) ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
  }
  if((calfH<=0)&&(calfL<=0))
  {
    snxt = kInfinity ;    // Travel perpendicular to z axis
    side = kNull;
  }
  // Radial Intersections
  //
  // Find intersection with cylinders at rmax/rmin
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //
  //            t1                t2                    t3

  t1   = 1.0 - v.z()*v.z() ;      // since v normalised
  t2   = p.x()*v.x() + p.y()*v.y() ;
  t3   = p.x()*p.x() + p.y()*p.y() ;

  if ( snxt > 10*(fDz+fRMax) )  { roi2 = 2*fRMax*fRMax; }
  else  { roi2 = snxt*snxt*t1 + 2*snxt*t2 + t3; }        // radius^2 on +-fDz

  if ( t1 > 0 ) // Check not parallel
  {
    // Calculate srd, r exit distance
     
    if ( (t2 >= 0.0) && (roi2 > fRMax*(fRMax + kRadTolerance)) )
    {
      // Delta r not negative => leaving via rmax

      deltaR = t3 - fRMax*fRMax ;

      // NOTE: Should use rho-fRMax<-kRadTolerance*0.5
      // - avoid sqrt for efficiency

      if ( deltaR < -kRadTolerance*fRMax )
      {
        b     = t2/t1 ;
        c     = deltaR/t1 ;
        d2    = b*b-c;
        if( d2 >= 0 ) { srd = c/( -b - std::sqrt(d2)); }
        else          { srd = 0.; }
        sider = kRMax ;
      }
      else
      {
        // On tolerant boundary & heading outwards (or perpendicular to)
        // outer radial surface -> leaving immediately

        if ( calcNorm ) 
        {
          *n         = G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0) ;
          *validNorm = true ;
        }
        return snxt = 0 ; // Leaving by rmax immediately
      }
    }             
    else if ( t2 < 0. ) // i.e.  t2 < 0; Possible rmin intersection
    {
      roMin2 = t3 - t2*t2/t1 ; // min ro2 of the plane of movement 

      if ( fRMin && (roMin2 < fRMin*(fRMin - kRadTolerance)) )
      {
        deltaR = t3 - fRMin*fRMin ;
        b      = t2/t1 ;
        c      = deltaR/t1 ;
        d2     = b*b - c ;

        if ( d2 >= 0 )   // Leaving via rmin
        {
          // NOTE: SHould use rho-rmin>kRadTolerance*0.5
          // - avoid sqrt for efficiency

          if (deltaR > kRadTolerance*fRMin)
          {
            srd = c/(-b+std::sqrt(d2)); 
            sider = kRMin ;
          }
          else
          {
            if ( calcNorm ) { *validNorm = false; }  // Concave side
            return snxt = 0.0;
          }
        }
        else    // No rmin intersect -> must be rmax intersect
        {
          deltaR = t3 - fRMax*fRMax ;
          c     = deltaR/t1 ;
          d2    = b*b-c;
          if( d2 >=0. )
          {
            srd    = -b + std::sqrt(d2) ;
            sider  = kRMax ;
          }
          else // Case: On the border+t2<kRadTolerance
               //       (v is perpendicular to the surface)
          {
            if (calcNorm)
            {
              *n = G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0) ;
              *validNorm = true ;
            }
            return snxt = 0.0;
          }
        }
      }
      else if ( roi2 > fRMax*(fRMax + kRadTolerance) )
           // No rmin intersect -> must be rmax intersect
      {
        deltaR = t3 - fRMax*fRMax ;
        b      = t2/t1 ;
        c      = deltaR/t1;
        d2     = b*b-c;
        if( d2 >= 0 )
        {
          srd    = -b + std::sqrt(d2) ;
          sider  = kRMax ;
        }
        else // Case: On the border+t2<kRadTolerance
             //       (v is perpendicular to the surface)
        {
          if (calcNorm)
          {
            *n = G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0) ;
            *validNorm = true ;
          }
          return snxt = 0.0;
        }
      }
    }
    // Phi Intersection

    if ( !fPhiFullCutTube )
    {
      // add angle calculation with correction 
      // of the difference in domain of atan2 and Sphi
      //
      vphi = std::atan2(v.y(),v.x()) ;
     
      if ( vphi < fSPhi - halfAngTolerance  )             { vphi += twopi; }
      else if ( vphi > fSPhi + fDPhi + halfAngTolerance ) { vphi -= twopi; }


      if ( p.x() || p.y() )  // Check if on z axis (rho not needed later)
      {
        // pDist -ve when inside

        pDistS = p.x()*sinSPhi - p.y()*cosSPhi ;
        pDistE = -p.x()*sinEPhi + p.y()*cosEPhi ;

        // Comp -ve when in direction of outwards normal

        compS   = -sinSPhi*v.x() + cosSPhi*v.y() ;
        compE   =  sinEPhi*v.x() - cosEPhi*v.y() ;
       
        sidephi = kNull;
        
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
              if( (std::fabs(xi)<=kCarTolerance)
               && (std::fabs(yi)<=kCarTolerance) )
              {
                sidephi = kSPhi;
                if (((fSPhi-halfAngTolerance)<=vphi)
                   &&((fSPhi+fDPhi+halfAngTolerance)>=vphi))
                {
                  sphi = kInfinity;
                }
              }
              else if ( yi*cosCPhi-xi*sinCPhi >=0 )
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
              
              if ((std::fabs(xi)<=kCarTolerance)&&(std::fabs(yi)<=kCarTolerance))
              {
                // Leaving via ending phi
                //
                if( !((fSPhi-halfAngTolerance <= vphi)
                     &&(fSPhi+fDPhi+halfAngTolerance >= vphi)) )
                {
                  sidephi = kEPhi ;
                  if ( pDistE <= -halfCarTolerance )  { sphi = sphi2 ; }
                  else                                { sphi = 0.0 ;   }
                }
              } 
              else    // Check intersecting with correct half-plane 

              if ( (yi*cosCPhi-xi*sinCPhi) >= 0)
              {
                // Leaving via ending phi
                //
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance ) { sphi = sphi2 ; }
                else                               { sphi = 0.0 ;   }
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
               
        if ( (fSPhi - halfAngTolerance <= vphi)
           && (vphi <= fSPhi + fDPhi + halfAngTolerance ) )
        {
          sphi = kInfinity ;
        }
        else
        {
          sidephi = kSPhi ; // arbitrary 
          sphi    = 0.0 ;
        }
      }
      if (sphi < snxt)  // Order intersecttions
      {
        snxt = sphi ;
        side = sidephi ;
      }
    }
    if (srd < snxt)  // Order intersections
    {
      snxt = srd ;
      side = sider ;
    }
  }
  if (calcNorm)
  {
    switch(side)
    {
      case kRMax:
        // Note: returned vector not normalised
        // (divide by fRMax for unit vector)
        //
        xi = p.x() + snxt*v.x() ;
        yi = p.y() + snxt*v.y() ;
        *n = G4ThreeVector(xi/fRMax,yi/fRMax,0) ;
        *validNorm = true ;
        break ;

      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;

      case kSPhi:
        if ( fDPhi <= pi )
        {
          *n         = G4ThreeVector(sinSPhi,-cosSPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kEPhi:
        if (fDPhi <= pi)
        {
          *n = G4ThreeVector(-sinEPhi,cosEPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kPZ:
        *n         = fHighNorm ;
        *validNorm = true ;
        break ;

      case kMZ:
        *n         = fLowNorm ;
        *validNorm = true ;
        break ;

      default:
        G4cout << G4endl ;
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
                << "snxt = "    << snxt/mm << " mm" << G4endl ;
        message.precision(oldprc) ;
        G4Exception("G4CutTubs::DistanceToOut(p,v,..)", "GeomSolids1002",
                    JustWarning, message);
        break ;
    }
  }
  if ( snxt<halfCarTolerance )  { snxt=0 ; }
  return snxt ;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4CutTubs::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safRMin,safRMax,safZLow,safZHigh,safePhi,safe,rho;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;  // Distance to R

  safRMin =  rho - fRMin ;
  safRMax =  fRMax - rho ;

  // Distances to ZCut(Low/High)

  // Dist to Low Cut
  //
  safZLow = std::fabs((p+vZ).dot(fLowNorm));

  // Dist to High Cut
  //
  safZHigh = std::fabs((p-vZ).dot(fHighNorm));
  safe = std::min(safZLow,safZHigh);

  if ( safRMin < safe ) { safe = safRMin; }
  if ( safRMax< safe )  { safe = safRMax; }

  // Check if phi divided, Calc distances closest phi plane
  //
  if ( !fPhiFullCutTube )
  {
    if ( p.y()*cosCPhi-p.x()*sinCPhi <= 0 )
    {
      safePhi = -(p.x()*sinSPhi - p.y()*cosSPhi) ;
    }
    else
    {
      safePhi = (p.x()*sinEPhi - p.y()*cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi ; }
  }
  if ( safe < 0 )  { safe = 0; }

  return safe ;
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

G4GeometryType G4CutTubs::GetEntityType() const
{
  return G4String("G4CutTubs");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4CutTubs::Clone() const
{
  return new G4CutTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4CutTubs::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4CutTubs\n"
     << " Parameters: \n"
     << "    inner radius : " << fRMin/mm << " mm \n"
     << "    outer radius : " << fRMax/mm << " mm \n"
     << "    half length Z: " << fDz/mm << " mm \n"
     << "    starting phi : " << fSPhi/degree << " degrees \n"
     << "    delta phi    : " << fDPhi/degree << " degrees \n"
     << "    low Norm     : " << fLowNorm     << "  \n" 
     << "    high Norm    : "  <<fHighNorm    << "  \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4CutTubs::GetPointOnSurface() const
{
  G4double xRand, yRand, zRand, phi, cosphi, sinphi, chose,
           aOne, aTwo, aThr, aFou;
  G4double rRand;
 
  aOne = 2.*fDz*fDPhi*fRMax;
  aTwo = 2.*fDz*fDPhi*fRMin;
  aThr = 0.5*fDPhi*(fRMax*fRMax-fRMin*fRMin);
  aFou = 2.*fDz*(fRMax-fRMin);

  phi    = G4RandFlat::shoot(fSPhi, fSPhi+fDPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);

  rRand  = GetRadiusInRing(fRMin,fRMax);
  
  if( (fSPhi == 0) && (fDPhi == twopi) ) { aFou = 0; }
  
  chose  = G4RandFlat::shoot(0.,aOne+aTwo+2.*aThr+2.*aFou);

  if( (chose >=0) && (chose < aOne) )
  {
    xRand = fRMax*cosphi;
    yRand = fRMax*sinphi;
    zRand = G4RandFlat::shoot(GetCutZ(G4ThreeVector(xRand,yRand,-fDz)),
                              GetCutZ(G4ThreeVector(xRand,yRand,fDz)));
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne) && (chose < aOne + aTwo) )
  {
    xRand = fRMin*cosphi;
    yRand = fRMin*sinphi;
    zRand = G4RandFlat::shoot(GetCutZ(G4ThreeVector(xRand,yRand,-fDz)),
                              GetCutZ(G4ThreeVector(xRand,yRand,fDz)));
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne + aTwo) && (chose < aOne + aTwo + aThr) )
  {
    xRand = rRand*cosphi;
    yRand = rRand*sinphi;
    zRand = GetCutZ(G4ThreeVector(xRand,yRand,fDz));
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne + aTwo + aThr) && (chose < aOne + aTwo + 2.*aThr) )
  {
    xRand = rRand*cosphi;
    yRand = rRand*sinphi;
    zRand = GetCutZ(G4ThreeVector(xRand,yRand,-fDz));
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne + aTwo + 2.*aThr)
        && (chose < aOne + aTwo + 2.*aThr + aFou) )
  {
    xRand = rRand*std::cos(fSPhi);
    yRand = rRand*std::sin(fSPhi);
    zRand = G4RandFlat::shoot(GetCutZ(G4ThreeVector(xRand,yRand,-fDz)),
                              GetCutZ(G4ThreeVector(xRand,yRand,fDz)));
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else
  {
    xRand = rRand*std::cos(fSPhi+fDPhi);
    yRand = rRand*std::sin(fSPhi+fDPhi);
    zRand = G4RandFlat::shoot(GetCutZ(G4ThreeVector(xRand,yRand,-fDz)),
                              GetCutZ(G4ThreeVector(xRand,yRand,fDz)));
    return G4ThreeVector  (xRand, yRand, zRand);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4CutTubs::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this) ;
}

G4Polyhedron* G4CutTubs::CreatePolyhedron () const 
{
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];

  G4Polyhedron *ph  = new G4Polyhedron;
  G4Polyhedron *ph1 = G4OTubs::CreatePolyhedron();
  G4int nn=ph1->GetNoVertices();
  G4int nf=ph1->GetNoFacets();
  G4double3* xyz = new G4double3[nn];  // number of nodes 
  G4int4*  faces = new G4int4[nf] ;    // number of faces

  for(G4int i=0;i<nn;i++)
  {
    xyz[i][0]=ph1->GetVertex(i+1).x();
    xyz[i][1]=ph1->GetVertex(i+1).y();
    G4double tmpZ=ph1->GetVertex(i+1).z();
    if(tmpZ>=fDz-kCarTolerance)
    {
      xyz[i][2]=GetCutZ(G4ThreeVector(xyz[i][0],xyz[i][1],fDz));
    }
    else if(tmpZ<=-fDz+kCarTolerance)
    {
      xyz[i][2]=GetCutZ(G4ThreeVector(xyz[i][0],xyz[i][1],-fDz));
    }
    else
    {
      xyz[i][2]=tmpZ;
    }
  }
  G4int iNodes[4];
  G4int *iEdge=0;
  G4int n;
  for(G4int i=0;i<nf;i++)
  {
    ph1->GetFacet(i+1,n,iNodes,iEdge);
    for(G4int k=0;k<n;k++)
    {
      faces[i][k]=iNodes[k];
    }
    for(G4int k=n;k<4;k++)
    {
      faces[i][k]=0;
    }
  }
  ph->createPolyhedron(nn,nf,xyz,faces);

  delete [] xyz;
  delete [] faces;
  delete ph1;

  return ph;
}

// Auxilary Methods for Solid
 
///////////////////////////////////////////////////////////////////////////
// Return true if Cutted planes are crossing 
// Check Intersection Points on OX and OY axes

G4bool G4CutTubs::IsCrossingCutPlanes() const
{
  G4double zXLow1,zXLow2,zYLow1,zYLow2;
  G4double zXHigh1,zXHigh2,zYHigh1,zYHigh2;

  zXLow1  = GetCutZ(G4ThreeVector(-fRMax,     0,-fDz));
  zXLow2  = GetCutZ(G4ThreeVector( fRMax,     0,-fDz));
  zYLow1  = GetCutZ(G4ThreeVector(     0,-fRMax,-fDz));
  zYLow2  = GetCutZ(G4ThreeVector(     0, fRMax,-fDz));
  zXHigh1 = GetCutZ(G4ThreeVector(-fRMax,     0, fDz));
  zXHigh2 = GetCutZ(G4ThreeVector( fRMax,     0, fDz));
  zYHigh1 = GetCutZ(G4ThreeVector(     0,-fRMax, fDz));
  zYHigh2 = GetCutZ(G4ThreeVector(     0, fRMax, fDz));
  if ( (zXLow1>zXHigh1) ||(zXLow2>zXHigh2)
    || (zYLow1>zYHigh1) ||(zYLow2>zYHigh2))  { return true; }

  return false;
}

///////////////////////////////////////////////////////////////////////////
//
// Return real Z coordinate of point on Cutted +/- fDZ plane

G4double G4CutTubs::GetCutZ(const G4ThreeVector& p) const
{
  G4double newz = p.z();  // p.z() should be either +fDz or -fDz
  if (p.z()<0)
  {
    if(fLowNorm.z()!=0.)
    {
       newz = -fDz-(p.x()*fLowNorm.x()+p.y()*fLowNorm.y())/fLowNorm.z();
    }
  }
  else
  {
    if(fHighNorm.z()!=0.)
    {
       newz = fDz-(p.x()*fHighNorm.x()+p.y()*fHighNorm.y())/fHighNorm.z();
    }
  }
  return newz;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate Min and Max Z for CutZ

void G4CutTubs::GetMaxMinZ(G4double& zmin,G4double& zmax)const

{
  G4double phiLow = std::atan2(fLowNorm.y(),fLowNorm.x());
  G4double phiHigh= std::atan2(fHighNorm.y(),fHighNorm.x());

  G4double xc=0, yc=0,z1;
  G4double z[8];
  G4bool in_range_low = false;
  G4bool in_range_hi = false;
 
  G4int i;
  for (i=0; i<2; i++)
  {
    if (phiLow<0)  { phiLow+=twopi; }
    G4double ddp = phiLow-fSPhi;
    if (ddp<0)  { ddp += twopi; }
    if (ddp <= fDPhi)
    {
      xc = fRMin*std::cos(phiLow);
      yc = fRMin*std::sin(phiLow);
      z1 = GetCutZ(G4ThreeVector(xc, yc, -fDz));
      xc = fRMax*std::cos(phiLow);
      yc = fRMax*std::sin(phiLow);
      z1 = std::min(z1, GetCutZ(G4ThreeVector(xc, yc, -fDz)));
      if (in_range_low)  { zmin = std::min(zmin, z1); }
      else               { zmin = z1; }
      in_range_low = true;
    }
    phiLow += pi;
    if (phiLow>twopi)  { phiLow-=twopi; }
  }
  for (i=0; i<2; i++)
  {
    if (phiHigh<0)  { phiHigh+=twopi; }
    G4double ddp = phiHigh-fSPhi;
    if (ddp<0)  { ddp += twopi; }
    if (ddp <= fDPhi)
    {
      xc = fRMin*std::cos(phiHigh);
      yc = fRMin*std::sin(phiHigh);
      z1 = GetCutZ(G4ThreeVector(xc, yc, fDz));
      xc = fRMax*std::cos(phiHigh);
      yc = fRMax*std::sin(phiHigh);
      z1 = std::min(z1, GetCutZ(G4ThreeVector(xc, yc, fDz)));
      if (in_range_hi)  { zmax = std::min(zmax, z1); }
      else              { zmax = z1; }
      in_range_hi = true;
    }
    phiHigh += pi;
    if (phiHigh>twopi)  { phiHigh-=twopi; }
  }

  xc = fRMin*std::cos(fSPhi);
  yc = fRMin*std::sin(fSPhi);
  z[0] = GetCutZ(G4ThreeVector(xc, yc, -fDz));
  z[4] = GetCutZ(G4ThreeVector(xc, yc, fDz));
 
  xc = fRMin*std::cos(fSPhi+fDPhi);
  yc = fRMin*std::sin(fSPhi+fDPhi);
  z[1] = GetCutZ(G4ThreeVector(xc, yc, -fDz));
  z[5] = GetCutZ(G4ThreeVector(xc, yc, fDz));
 
  xc = fRMax*std::cos(fSPhi);
  yc = fRMax*std::sin(fSPhi);
  z[2] = GetCutZ(G4ThreeVector(xc, yc, -fDz));
  z[6] = GetCutZ(G4ThreeVector(xc, yc, fDz));
 
  xc = fRMax*std::cos(fSPhi+fDPhi);
  yc = fRMax*std::sin(fSPhi+fDPhi);
  z[3] = GetCutZ(G4ThreeVector(xc, yc, -fDz));
  z[7] = GetCutZ(G4ThreeVector(xc, yc, fDz));
 
  // Find min/max

  z1=z[0];
  for (i = 1; i < 4; i++)
  {
    if(z[i] < z[i-1])z1=z[i];
  }
    
  if (in_range_low)
  {
    zmin = std::min(zmin, z1);
  }
  else
  {
    zmin = z1;
  }
  z1=z[4];
  for (i = 1; i < 4; i++)
  {
    if(z[4+i] > z[4+i-1])  { z1=z[4+i]; }
  }

  if (in_range_hi)  { zmax = std::max(zmax, z1); }
  else              { zmax = z1; }
}
