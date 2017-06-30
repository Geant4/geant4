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
// class G4OTubs
//
//   Temporary copy of original G4Tubs code for use by G4CutTubs.
//
/////////////////////////////////////////////////////////////////////////

#include "G4OTubs.hh"

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4OTubs::G4OTubs( const G4String &pName,
                      G4double pRMin, G4double pRMax,
                      G4double pDz,
                      G4double pSPhi, G4double pDPhi )
  : G4CSGSolid(pName), fRMin(pRMin), fRMax(pRMax), fDz(pDz), fSPhi(0), fDPhi(0)
{

  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  halfCarTolerance=kCarTolerance*0.5;
  halfRadTolerance=kRadTolerance*0.5;
  halfAngTolerance=kAngTolerance*0.5;

  if (pDz<=0) // Check z-len
  {
    std::ostringstream message;
    message << "Negative Z half-length (" << pDz << ") in solid: " << GetName();
    G4Exception("G4Tubs::G4Tubs()", "GeomSolids0002", FatalException, message);
  }
  if ( (pRMin >= pRMax) || (pRMin < 0) ) // Check radii
  {
    std::ostringstream message;
    message << "Invalid values for radii in solid: " << GetName()
            << G4endl
            << "        pRMin = " << pRMin << ", pRMax = " << pRMax;
    G4Exception("G4Tubs::G4Tubs()", "GeomSolids0002", FatalException, message);
  }

  // Check angles
  //
  CheckPhiAngles(pSPhi, pDPhi);
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4OTubs::G4OTubs( __void__& a )
  : G4CSGSolid(a), kRadTolerance(0.), kAngTolerance(0.),
    fRMin(0.), fRMax(0.), fDz(0.), fSPhi(0.), fDPhi(0.),
    sinCPhi(0.), cosCPhi(0.), cosHDPhiOT(0.), cosHDPhiIT(0.),
    sinSPhi(0.), cosSPhi(0.), sinEPhi(0.), cosEPhi(0.),
    fPhiFullTube(false), halfCarTolerance(0.), halfRadTolerance(0.),
    halfAngTolerance(0.)

{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4OTubs::~G4OTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4OTubs::G4OTubs(const G4OTubs& rhs)
  : G4CSGSolid(rhs),
    kRadTolerance(rhs.kRadTolerance), kAngTolerance(rhs.kAngTolerance),
    fRMin(rhs.fRMin), fRMax(rhs.fRMax), fDz(rhs.fDz),
    fSPhi(rhs.fSPhi), fDPhi(rhs.fDPhi),
    sinCPhi(rhs.sinCPhi), cosCPhi(rhs.cosCPhi),
    cosHDPhiOT(rhs.cosHDPhiOT), cosHDPhiIT(rhs.cosHDPhiIT),
    sinSPhi(rhs.sinSPhi), cosSPhi(rhs.cosSPhi),
    sinEPhi(rhs.sinEPhi), cosEPhi(rhs.cosEPhi), fPhiFullTube(rhs.fPhiFullTube),
    halfCarTolerance(rhs.halfCarTolerance),
    halfRadTolerance(rhs.halfRadTolerance),
    halfAngTolerance(rhs.halfAngTolerance)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4OTubs& G4OTubs::operator = (const G4OTubs& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   kRadTolerance = rhs.kRadTolerance; kAngTolerance = rhs.kAngTolerance;
   fRMin = rhs.fRMin; fRMax = rhs.fRMax; fDz = rhs.fDz;
   fSPhi = rhs.fSPhi; fDPhi = rhs.fDPhi;
   sinCPhi = rhs.sinCPhi; cosCPhi = rhs.cosCPhi;
   cosHDPhiOT = rhs.cosHDPhiOT; cosHDPhiIT = rhs.cosHDPhiIT;
   sinSPhi = rhs.sinSPhi; cosSPhi = rhs.cosSPhi;
   sinEPhi = rhs.sinEPhi; cosEPhi = rhs.cosEPhi;
   fPhiFullTube = rhs.fPhiFullTube;
   halfCarTolerance = rhs.halfCarTolerance;
   halfRadTolerance = rhs.halfRadTolerance;
   halfAngTolerance = rhs.halfAngTolerance;

   return *this;
}

/////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4OTubs::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
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
    G4Exception("G4OTubs::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4OTubs::CalculateExtent( const EAxis              pAxis,
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

///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4OTubs::Inside( const G4ThreeVector& p ) const
{
  G4double r2,pPhi,tolRMin,tolRMax;
  EInside in = kOutside ;

  if (std::fabs(p.z()) <= fDz - halfCarTolerance)
  {
    r2 = p.x()*p.x() + p.y()*p.y() ;

    if (fRMin) { tolRMin = fRMin + halfRadTolerance ; }
    else       { tolRMin = 0 ; }

    tolRMax = fRMax - halfRadTolerance ;
      
    if ((r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax))
    {
      if ( fPhiFullTube )
      {
        in = kInside ;
      }
      else
      {
        // Try inner tolerant phi boundaries (=>inside)
        // if not inside, try outer tolerant phi boundaries

        if ( (tolRMin==0) && (std::fabs(p.x())<=halfCarTolerance)
                          && (std::fabs(p.y())<=halfCarTolerance) )
        {
          in=kSurface;
        }
        else
        {
          pPhi = std::atan2(p.y(),p.x()) ;
          if ( pPhi < -halfAngTolerance )  { pPhi += twopi; } // 0<=pPhi<2pi

          if ( fSPhi >= 0 )
          {
            if ( (std::fabs(pPhi) < halfAngTolerance)
              && (std::fabs(fSPhi + fDPhi - twopi) < halfAngTolerance) )
            { 
              pPhi += twopi ; // 0 <= pPhi < 2pi
            }
            if ( (pPhi >= fSPhi + halfAngTolerance)
              && (pPhi <= fSPhi + fDPhi - halfAngTolerance) )
            {
              in = kInside ;
            }
            else if ( (pPhi >= fSPhi - halfAngTolerance)
                   && (pPhi <= fSPhi + fDPhi + halfAngTolerance) )
            {
              in = kSurface ;
            }
          }
          else  // fSPhi < 0
          {
            if ( (pPhi <= fSPhi + twopi - halfAngTolerance)
              && (pPhi >= fSPhi + fDPhi  + halfAngTolerance) ) {;} //kOutside
            else if ( (pPhi <= fSPhi + twopi + halfAngTolerance)
                   && (pPhi >= fSPhi + fDPhi  - halfAngTolerance) )
            {
              in = kSurface ;
            }
            else
            {
              in = kInside ;
            }
          }
        }                    
      }
    }
    else  // Try generous boundaries
    {
      tolRMin = fRMin - halfRadTolerance ;
      tolRMax = fRMax + halfRadTolerance ;

      if ( tolRMin < 0 )  { tolRMin = 0; }

      if ( (r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax) )
      {
        if (fPhiFullTube || (r2 <=halfRadTolerance*halfRadTolerance) )
        {                        // Continuous in phi or on z-axis
          in = kSurface ;
        }
        else // Try outer tolerant phi boundaries only
        {
          pPhi = std::atan2(p.y(),p.x()) ;

          if ( pPhi < -halfAngTolerance)  { pPhi += twopi; } // 0<=pPhi<2pi
          if ( fSPhi >= 0 )
          {
            if ( (std::fabs(pPhi) < halfAngTolerance)
              && (std::fabs(fSPhi + fDPhi - twopi) < halfAngTolerance) )
            { 
              pPhi += twopi ; // 0 <= pPhi < 2pi
            }
            if ( (pPhi >= fSPhi - halfAngTolerance)
              && (pPhi <= fSPhi + fDPhi + halfAngTolerance) )
            {
              in = kSurface ;
            }
          }
          else  // fSPhi < 0
          {
            if ( (pPhi <= fSPhi + twopi - halfAngTolerance)
              && (pPhi >= fSPhi + fDPhi + halfAngTolerance) ) {;} // kOutside
            else
            {
              in = kSurface ;
            }
          }
        }
      }
    }
  }
  else if (std::fabs(p.z()) <= fDz + halfCarTolerance)
  {                                          // Check within tolerant r limits
    r2      = p.x()*p.x() + p.y()*p.y() ;
    tolRMin = fRMin - halfRadTolerance ;
    tolRMax = fRMax + halfRadTolerance ;

    if ( tolRMin < 0 )  { tolRMin = 0; }

    if ( (r2 >= tolRMin*tolRMin) && (r2 <= tolRMax*tolRMax) )
    {
      if (fPhiFullTube || (r2 <=halfRadTolerance*halfRadTolerance))
      {                        // Continuous in phi or on z-axis
        in = kSurface ;
      }
      else // Try outer tolerant phi boundaries
      {
        pPhi = std::atan2(p.y(),p.x()) ;

        if ( pPhi < -halfAngTolerance )  { pPhi += twopi; }  // 0<=pPhi<2pi
        if ( fSPhi >= 0 )
        {
          if ( (std::fabs(pPhi) < halfAngTolerance)
            && (std::fabs(fSPhi + fDPhi - twopi) < halfAngTolerance) )
          { 
            pPhi += twopi ; // 0 <= pPhi < 2pi
          }
          if ( (pPhi >= fSPhi - halfAngTolerance)
            && (pPhi <= fSPhi + fDPhi + halfAngTolerance) )
          {
            in = kSurface;
          }
        }
        else  // fSPhi < 0
        {
          if ( (pPhi <= fSPhi + twopi - halfAngTolerance)
            && (pPhi >= fSPhi + fDPhi  + halfAngTolerance) ) {;}
          else
          {
            in = kSurface ;
          }
        }      
      }
    }
  }
  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4OTubs::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4int noSurfaces = 0;
  G4double rho, pPhi;
  G4double distZ, distRMin, distRMax;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;

  G4ThreeVector norm, sumnorm(0.,0.,0.);
  G4ThreeVector nZ = G4ThreeVector(0, 0, 1.0);
  G4ThreeVector nR, nPs, nPe;

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y());

  distRMin = std::fabs(rho - fRMin);
  distRMax = std::fabs(rho - fRMax);
  distZ    = std::fabs(std::fabs(p.z()) - fDz);

  if (!fPhiFullTube)    // Protected against (0,0,z) 
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
  if (distZ <= halfCarTolerance)  
  {
    noSurfaces ++;
    if ( p.z() >= 0.)  { sumnorm += nZ; }
    else               { sumnorm -= nZ; }
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4Tubs::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" );
    G4int oldprc = G4cout.precision(20);
    G4cout<< "G4Tubs::SN ( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "
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

G4ThreeVector G4OTubs::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  ENorm side ;
  G4ThreeVector norm ;
  G4double rho, phi ;
  G4double distZ, distRMin, distRMax, distSPhi, distEPhi, distMin ;

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

  distRMin = std::fabs(rho - fRMin) ;
  distRMax = std::fabs(rho - fRMax) ;
  distZ    = std::fabs(std::fabs(p.z()) - fDz) ;

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
  if (!fPhiFullTube  &&  rho ) // Protected against (0,0,z) 
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
      if ( p.z() > 0 )  { norm = G4ThreeVector(0,0,1) ; }
      else              { norm = G4ThreeVector(0,0,-1); }
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
      G4Exception("G4Tubs::ApproxSurfaceNormal()",
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

G4double G4OTubs::DistanceToIn( const G4ThreeVector& p,
                               const G4ThreeVector& v  ) const
{
  G4double snxt = kInfinity ;      // snxt = default return value
  G4double tolORMin2, tolIRMax2 ;  // 'generous' radii squared
  G4double tolORMax2, tolIRMin2, tolODz, tolIDz ;
  const G4double dRmax = 100.*fRMax;

  // Intersection point variables
  //
  G4double Dist, sd, xi, yi, zi, rho2, inum, iden, cosPsi, Comp ;
  G4double t1, t2, t3, b, c, d ;     // Quadratic solver variables 
  
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

  // Intersection with Z surfaces

  tolIDz = fDz - halfCarTolerance ;
  tolODz = fDz + halfCarTolerance ;

  if (std::fabs(p.z()) >= tolIDz)
  {
    if ( p.z()*v.z() < 0 )    // at +Z going in -Z or visa versa
    {
      sd = (std::fabs(p.z()) - fDz)/std::fabs(v.z()) ;  // Z intersect distance

      if(sd < 0.0)  { sd = 0.0; }

      xi   = p.x() + sd*v.x() ;                // Intersection coords
      yi   = p.y() + sd*v.y() ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
        if (!fPhiFullTube && rho2)
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
      if ( snxt<halfCarTolerance )  { snxt=0; }
      return snxt ;  // On/outside extent, and heading away
                     // -> cannot intersect
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
      // Try outer cylinder intersection
      //          c=(t3-fRMax*fRMax)/t1;

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
          if (std::fabs(zi)<=tolODz)
          {
            // Z ok. Check phi intersection if reqd
            //
            if (fPhiFullTube)
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
        }    //  end if (sd>=0)
      }      //  end if (d>=0)
    }        //  end if (r>=fRMax)
    else 
    {
      // Inside outer radius :
      // check not inside, and heading through tubs (-> 0 to in)

      if ((t3 > tolIRMin2) && (t2 < 0) && (std::fabs(p.z()) <= tolIDz))
      {
        // Inside both radii, delta r -ve, inside z extent

        if (!fPhiFullTube)
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
        } // end if   (!fPhiFullTube)
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
        if (sd >= -halfCarTolerance)  // check forwards
        {
          // Check z intersection
          //
          if(sd < 0.0)  { sd = 0.0; }
          if ( sd>dRmax ) // Avoid rounding errors due to precision issues seen
          {               // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd-std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
          zi = p.z() + sd*v.z() ;
          if (std::fabs(zi) <= tolODz)
          {
            // Z ok. Check phi
            //
            if ( fPhiFullTube )
            {
              return sd ; 
            }
            else
            {
              xi     = p.x() + sd*v.x() ;
              yi     = p.y() + sd*v.y() ;
              cosPsi = (xi*cosCPhi + yi*sinCPhi)/fRMin ;
              if (cosPsi >= cosHDPhiIT)
              {
                // Good inner radius isect
                // - but earlier phi isect still possible

                snxt = sd ;
              }
            }
          }        //    end if std::fabs(zi)
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
  if ( !fPhiFullTube )
  {
    // First phi surface (Starting phi)
    //
    Comp    = v.x()*sinSPhi - v.y()*cosSPhi ;
                    
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
          if ( std::fabs(zi) <= tolODz )
          {
            xi   = p.x() + sd*v.x() ;
            yi   = p.y() + sd*v.y() ;
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
          }
        }
      }    
    }
      
    // Second phi surface (Ending phi)

    Comp    = -(v.x()*sinEPhi - v.y()*cosEPhi) ;
        
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
          if ( std::fabs(zi) <= tolODz )
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
              if ( (yi*cosCPhi-xi*sinCPhi) >= 0 ) { snxt = sd; }
            }                         //?? >=-halfCarTolerance
          }
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

G4double G4OTubs::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe=0.0, rho, safe1, safe2, safe3 ;
  G4double safePhi, cosPsi ;

  rho   = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;
  safe1 = fRMin - rho ;
  safe2 = rho - fRMax ;
  safe3 = std::fabs(p.z()) - fDz ;

  if ( safe1 > safe2 ) { safe = safe1; }
  else                 { safe = safe2; }
  if ( safe3 > safe )  { safe = safe3; }

  if ( (!fPhiFullTube) && (rho) )
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

G4double G4OTubs::DistanceToOut( const G4ThreeVector& p,
                                const G4ThreeVector& v,
                                const G4bool calcNorm,
                                      G4bool *validNorm,
                                      G4ThreeVector *n    ) const
{  
  ESide side=kNull , sider=kNull, sidephi=kNull ;
  G4double snxt, srd=kInfinity, sphi=kInfinity, pdist ;
  G4double deltaR, t1, t2, t3, b, c, d2, roMin2 ;

  // Vars for phi intersection:

  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, vphi, roi2 ;
 
  // Z plane intersection

  if (v.z() > 0 )
  {
    pdist = fDz - p.z() ;
    if ( pdist > halfCarTolerance )
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
      return snxt = 0 ;
    }
  }
  else if ( v.z() < 0 )
  {
    pdist = fDz + p.z() ;

    if ( pdist > halfCarTolerance )
    {
      snxt = -pdist/v.z() ;
      side = kMZ ;
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
  else
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
            srd     = -b + std::sqrt(d2) ;
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
          srd     = -b + std::sqrt(d2) ;
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

    if ( !fPhiFullTube )
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
              
              if ( (std::fabs(xi)<=kCarTolerance)
                && (std::fabs(yi)<=kCarTolerance))
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
        G4Exception("G4Tubs::DistanceToOut(p,v,..)", "GeomSolids1002",
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

G4double G4OTubs::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0, rho, safeR1, safeR2, safeZ, safePhi ;
  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
    G4int oldprc = G4cout.precision(16) ;
    G4cout << G4endl ;
    DumpInfo();
    G4cout << "Position:"  << G4endl << G4endl ;
    G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
    G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
    G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
    G4cout.precision(oldprc) ;
    G4Exception("G4Tubs::DistanceToOut(p)", "GeomSolids1002",
                JustWarning, "Point p is outside !?");
  }
#endif

  if ( fRMin )
  {
    safeR1 = rho   - fRMin ;
    safeR2 = fRMax - rho ;
 
    if ( safeR1 < safeR2 ) { safe = safeR1 ; }
    else                   { safe = safeR2 ; }
  }
  else
  {
    safe = fRMax - rho ;
  }
  safeZ = fDz - std::fabs(p.z()) ;

  if ( safeZ < safe )  { safe = safeZ ; }

  // Check if phi divided, Calc distances closest phi plane
  //
  if ( !fPhiFullTube )
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
  if ( safe < 0 )  { safe = 0 ; }

  return safe ;  
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

G4GeometryType G4OTubs::GetEntityType() const
{
  return G4String("G4Tubs");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4OTubs::Clone() const
{
  return new G4OTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4OTubs::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Tubs\n"
     << " Parameters: \n"
     << "    inner radius : " << fRMin/mm << " mm \n"
     << "    outer radius : " << fRMax/mm << " mm \n"
     << "    half length Z: " << fDz/mm << " mm \n"
     << "    starting phi : " << fSPhi/degree << " degrees \n"
     << "    delta phi    : " << fDPhi/degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4OTubs::GetPointOnSurface() const
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
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne) && (chose < aOne + aTwo) )
  {
    xRand = fRMin*cosphi;
    yRand = fRMin*sinphi;
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne + aTwo) && (chose < aOne + aTwo + aThr) )
  {
    xRand = rRand*cosphi;
    yRand = rRand*sinphi;
    zRand = fDz;
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne + aTwo + aThr) && (chose < aOne + aTwo + 2.*aThr) )
  {
    xRand = rRand*cosphi;
    yRand = rRand*sinphi;
    zRand = -1.*fDz;
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else if( (chose >= aOne + aTwo + 2.*aThr)
        && (chose < aOne + aTwo + 2.*aThr + aFou) )
  {
    xRand = rRand*std::cos(fSPhi);
    yRand = rRand*std::sin(fSPhi);
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector  (xRand, yRand, zRand);
  }
  else
  {
    xRand = rRand*std::cos(fSPhi+fDPhi);
    yRand = rRand*std::sin(fSPhi+fDPhi);
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector  (xRand, yRand, zRand);
  }
}

///////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4OTubs::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this) ;
}

G4Polyhedron* G4OTubs::CreatePolyhedron () const 
{
  return new G4PolyhedronTubs (fRMin, fRMax, fDz, fSPhi, fDPhi) ;
}
