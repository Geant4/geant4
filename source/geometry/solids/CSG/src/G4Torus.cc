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
// $Id: G4Torus.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// class G4Torus
//
// Implementation
//
// 16.12.16 H.Burkhardt: use radius differences and hypot to improve precision
// 28.10.16 E.Tcherniaev: reimplemented CalculateExtent(),
//                      removed CreateRotatedVertices()
// 05.04.12 M.Kelsey:   Use sqrt(r) in GetPointOnSurface() for uniform points
// 02.10.07 T.Nikitina: Bug fixed in SolveNumericJT(), b.969:segmentation fault.
//                      rootsrefined is used only if the number of refined roots
//                      is the same as for primary roots. 
// 02.10.07 T.Nikitina: Bug fixed in CalculateExtent() for case of non-rotated
//                      full-phi torus:protect against negative value for sqrt,
//                      correct  formula for delta.  
// 20.11.05 V.Grichine: Bug fixed in Inside(p) for phi sections, b.810 
// 25.08.05 O.Link: new methods for DistanceToIn/Out using JTPolynomialSolver
// 07.06.05 V.Grichine: SurfaceNormal(p) for rho=0, Constructor as G4Cons 
// 03.05.05 V.Grichine: SurfaceNormal(p) according to J. Apostolakis proposal
// 18.03.04 V.Grichine: bug fixed in DistanceToIn(p)
// 11.01.01 E.Medernach: Use G4PolynomialSolver to find roots
// 03.10.00 E.Medernach: SafeNewton added
// 31.08.00 E.Medernach: numerical computation of roots wuth bounding
//                       volume technique
// 26.05.00 V.Grichine: new fuctions developed by O.Cremonesi were added
// 06.03.00 V.Grichine: modifications in Distance ToOut(p,v,...)
// 19.11.99 V.Grichine: side = kNull in Distance ToOut(p,v,...)
// 09.10.98 V.Grichine: modifications in Distance ToOut(p,v,...)
// 30.10.96 V.Grichine: first implementation with G4Tubs elements in Fs
//

#include "G4Torus.hh"

#if !(defined(G4GEOM_USE_UTORUS) && defined(G4GEOM_USE_SYS_USOLIDS))

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"
#include "G4JTPolynomialSolver.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

using namespace CLHEP;

///////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4Torus::G4Torus( const G4String &pName,
                        G4double pRmin,
                        G4double pRmax,
                        G4double pRtor,
                        G4double pSPhi,
                        G4double pDPhi)
  : G4CSGSolid(pName)
{
  SetAllParameters(pRmin, pRmax, pRtor, pSPhi, pDPhi);
}

////////////////////////////////////////////////////////////////////////////
//
//

void
G4Torus::SetAllParameters( G4double pRmin,
                           G4double pRmax,
                           G4double pRtor,
                           G4double pSPhi,
                           G4double pDPhi )
{
  const G4double fEpsilon = 4.e-11;  // relative tolerance of radii

  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fRebuildPolyhedron = true;

  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  halfCarTolerance = 0.5*kCarTolerance;
  halfAngTolerance = 0.5*kAngTolerance;

  if ( pRtor >= pRmax+1.e3*kCarTolerance )  // Check swept radius, as in G4Cons
  {
    fRtor = pRtor ;
  }
  else
  {
    std::ostringstream message;
    message << "Invalid swept radius for Solid: " << GetName() << G4endl
            << "        pRtor = " << pRtor << ", pRmax = " << pRmax;
    G4Exception("G4Torus::SetAllParameters()",
                "GeomSolids0002", FatalException, message);
  }

  // Check radii, as in G4Cons
  //
  if ( pRmin < pRmax - 1.e2*kCarTolerance && pRmin >= 0 )
  {
    if (pRmin >= 1.e2*kCarTolerance) { fRmin = pRmin ; }
    else                             { fRmin = 0.0   ; }
    fRmax = pRmax ;
  }
  else
  {
    std::ostringstream message;
    message << "Invalid values of radii for Solid: " << GetName() << G4endl
            << "        pRmin = " << pRmin << ", pRmax = " << pRmax;
    G4Exception("G4Torus::SetAllParameters()",
                "GeomSolids0002", FatalException, message);
  }

  // Relative tolerances
  //
  fRminTolerance = (fRmin)
                 ? 0.5*std::max( kRadTolerance, fEpsilon*(fRtor-fRmin )) : 0;
  fRmaxTolerance = 0.5*std::max( kRadTolerance, fEpsilon*(fRtor+fRmax) );

  // Check angles
  //
  if ( pDPhi >= twopi )  { fDPhi = twopi ; }
  else
  {
    if (pDPhi > 0)       { fDPhi = pDPhi ; }
    else
    {
      std::ostringstream message;
      message << "Invalid Z delta-Phi for Solid: " << GetName() << G4endl
              << "        pDPhi = " << pDPhi;
      G4Exception("G4Torus::SetAllParameters()",
                  "GeomSolids0002", FatalException, message);
    }
  }
  
  // Ensure psphi in 0-2PI or -2PI-0 range if shape crosses 0
  //
  fSPhi = pSPhi;

  if (fSPhi < 0)  { fSPhi = twopi-std::fmod(std::fabs(fSPhi),twopi) ; }
  else            { fSPhi = std::fmod(fSPhi,twopi) ; }

  if (fSPhi+fDPhi > twopi)  { fSPhi-=twopi ; }
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Torus::G4Torus( __void__& a )
  : G4CSGSolid(a), fRmin(0.), fRmax(0.), fRtor(0.), fSPhi(0.),
    fDPhi(0.), fRminTolerance(0.), fRmaxTolerance(0. ),
    kRadTolerance(0.), kAngTolerance(0.),
    halfCarTolerance(0.), halfAngTolerance(0.)
{
}

//////////////////////////////////////////////////////////////////////
//
// Destructor

G4Torus::~G4Torus()
{}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Torus::G4Torus(const G4Torus& rhs)
  : G4CSGSolid(rhs), fRmin(rhs.fRmin),fRmax(rhs.fRmax),
    fRtor(rhs.fRtor),fSPhi(rhs.fSPhi),fDPhi(rhs.fDPhi),
    fRminTolerance(rhs.fRminTolerance), fRmaxTolerance(rhs.fRmaxTolerance),
    kRadTolerance(rhs.kRadTolerance), kAngTolerance(rhs.kAngTolerance),
    halfCarTolerance(rhs.halfCarTolerance),
    halfAngTolerance(rhs.halfAngTolerance)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Torus& G4Torus::operator = (const G4Torus& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fRmin = rhs.fRmin; fRmax = rhs.fRmax;
   fRtor = rhs.fRtor; fSPhi = rhs.fSPhi; fDPhi = rhs.fDPhi;
   fRminTolerance = rhs.fRminTolerance; fRmaxTolerance = rhs.fRmaxTolerance;
   kRadTolerance = rhs.kRadTolerance; kAngTolerance = rhs.kAngTolerance;
   halfCarTolerance = rhs.halfCarTolerance;
   halfAngTolerance = rhs.halfAngTolerance;

   return *this;
}

//////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Torus::ComputeDimensions(       G4VPVParameterisation* p,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}



////////////////////////////////////////////////////////////////////////////////
//
// Calculate the real roots to torus surface. 
// Returns negative solutions as well.

void G4Torus::TorusRootsJT( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                                  G4double r,
                                  std::vector<G4double>& roots ) const
{

  G4int i, num ;
  G4double c[5], srd[4], si[4] ;

  G4double Rtor2 = fRtor*fRtor, r2 = r*r  ;

  G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
  G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;

  G4double d=pRad2 - Rtor2;
  c[0] = 1.0 ;
  c[1] = 4*pDotV ;
  c[2] = 2*( (d + 2*pDotV*pDotV  - r2) + 2*Rtor2*v.z()*v.z());
  c[3] = 4*(pDotV*(d - r2) + 2*Rtor2*p.z()*v.z()) ;
  c[4] = (d-r2)*(d-r2) +4*Rtor2*(p.z()*p.z()-r2);

  G4JTPolynomialSolver  torusEq;

  num = torusEq.FindRoots( c, 4, srd, si );
  
  for ( i = 0; i < num; i++ ) 
  {
    if( si[i] == 0. )  { roots.push_back(srd[i]) ; }  // store real roots
  }  

  std::sort(roots.begin() , roots.end() ) ;  // sorting  with <
}

//////////////////////////////////////////////////////////////////////////////
//
// Interface for DistanceToIn and DistanceToOut.
// Calls TorusRootsJT and returns the smalles possible distance to 
// the surface.
// Attention: Difference in DistanceToIn/Out for points p on the surface.

G4double G4Torus::SolveNumericJT( const G4ThreeVector& p,
                                  const G4ThreeVector& v,
                                        G4double r,
                                        G4bool IsDistanceToIn ) const
{
  G4double bigdist = 10*mm ;
  G4double tmin = kInfinity ;
  G4double t, scal ;

  // calculate the distances to the intersections with the Torus
  // from a given point p and direction v.
  //
  std::vector<G4double> roots ;
  std::vector<G4double> rootsrefined ;
  TorusRootsJT(p,v,r,roots) ;

  G4ThreeVector ptmp ;

  // determine the smallest non-negative solution
  //
  for ( size_t k = 0 ; k<roots.size() ; k++ )
  {
    t = roots[k] ;

    if ( t < -halfCarTolerance )  { continue ; }  // skip negative roots

    if ( t > bigdist && t<kInfinity )    // problem with big distances
    {
      ptmp = p + t*v ;
      TorusRootsJT(ptmp,v,r,rootsrefined) ;
      if ( rootsrefined.size()==roots.size() )
      {
        t = t + rootsrefined[k] ;
      }
    }

    ptmp = p + t*v ;   // calculate the position of the proposed intersection

    G4double theta = std::atan2(ptmp.y(),ptmp.x());
    
    if ( fSPhi >= 0 )
    {
      if ( theta < - halfAngTolerance )  { theta += twopi; }
      if ( (std::fabs(theta) < halfAngTolerance)
        && (std::fabs(fSPhi + fDPhi - twopi) < halfAngTolerance) )
      { 
        theta += twopi ; // 0 <= theta < 2pi
      }
    }
    if ((fSPhi <= -pi )&&(theta>halfAngTolerance)) { theta = theta-twopi; }
       
    // We have to verify if this root is inside the region between
    // fSPhi and fSPhi + fDPhi
    //
    if ( (theta - fSPhi >= - halfAngTolerance)
      && (theta - (fSPhi + fDPhi) <=  halfAngTolerance) )
    {
      // check if P is on the surface, and called from DistanceToIn
      // DistanceToIn has to return 0.0 if particle is going inside the solid

      if ( IsDistanceToIn == true )
      {
        if (std::fabs(t) < halfCarTolerance )
        {
          // compute scalar product at position p : v.n
          // ( n taken from SurfaceNormal, not normalized )

          scal = v* G4ThreeVector( p.x()*(1-fRtor/std::hypot(p.x(),p.y())),
                                   p.y()*(1-fRtor/std::hypot(p.x(),p.y())),
                                   p.z() );

          // change sign in case of inner radius
          //
          if ( r == GetRmin() )  { scal = -scal ; }
          if ( scal < 0 )  { return 0.0  ; }
        }
      }

      // check if P is on the surface, and called from DistanceToOut
      // DistanceToIn has to return 0.0 if particle is leaving the solid

      if ( IsDistanceToIn == false )
      {
        if (std::fabs(t) < halfCarTolerance )
        {
          // compute scalar product at position p : v.n   
          //
          scal = v* G4ThreeVector( p.x()*(1-fRtor/std::hypot(p.x(),p.y())),
                                   p.y()*(1-fRtor/std::hypot(p.x(),p.y())),
                                   p.z() );

          // change sign in case of inner radius
          //
          if ( r == GetRmin() )  { scal = -scal ; }
          if ( scal > 0 )  { return 0.0  ; }
        }
      }

      // check if distance is larger than 1/2 kCarTolerance
      //
      if(  t > halfCarTolerance  )
      {
        tmin = t  ;
        return tmin  ;
      }
    }
  }

  return tmin;
}

/////////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Torus::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
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
    G4Exception("G4Torus::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Torus::CalculateExtent( const EAxis pAxis,
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

//////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Torus::Inside( const G4ThreeVector& p ) const
{
  G4double r, pt2, pPhi, tolRMin, tolRMax ;

  EInside in = kOutside ;

  // General precals
  //
  r   = std::hypot(p.x(),p.y());
  pt2 = p.z()*p.z() + (r-fRtor)*(r-fRtor);

  if (fRmin) tolRMin = fRmin + fRminTolerance ;
  else       tolRMin = 0 ;

  tolRMax = fRmax - fRmaxTolerance;
      
  if (pt2 >= tolRMin*tolRMin && pt2 <= tolRMax*tolRMax )
  {
    if ( fDPhi == twopi || pt2 == 0 )  // on torus swept axis
    {
      in = kInside ;
    }
    else
    {
      // Try inner tolerant phi boundaries (=>inside)
      // if not inside, try outer tolerant phi boundaries

      pPhi = std::atan2(p.y(),p.x()) ;

      if ( pPhi < -halfAngTolerance )  { pPhi += twopi ; }  // 0<=pPhi<2pi
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
            && (pPhi >= fSPhi + fDPhi  + halfAngTolerance) )  {;}
          else
          {
            in = kSurface ;
          }
      }
    }
  }
  else   // Try generous boundaries
  {
    tolRMin = fRmin - fRminTolerance ;
    tolRMax = fRmax + fRmaxTolerance ;

    if (tolRMin < 0 )  { tolRMin = 0 ; }

    if ( (pt2 >= tolRMin*tolRMin) && (pt2 <= tolRMax*tolRMax) )
    {
      if ( (fDPhi == twopi) || (pt2 == 0) ) // Continuous in phi or on z-axis
      {
        in = kSurface ;
      }
      else // Try outer tolerant phi boundaries only
      {
        pPhi = std::atan2(p.y(),p.x()) ;

        if ( pPhi < -halfAngTolerance )  { pPhi += twopi ; }  // 0<=pPhi<2pi
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
            && (pPhi >= fSPhi + fDPhi  + halfAngTolerance) )  {;}
          else
          {
            in = kSurface ;
          }
        }
      }
    }
  }
  return in ;
}

/////////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Torus::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4int noSurfaces = 0;  
  G4double rho, pt, pPhi;
  G4double distRMin = kInfinity;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;

  // To cope with precision loss
  //
  const G4double delta = std::max(10.0*kCarTolerance,
                                  1.0e-8*(fRtor+fRmax));
  const G4double dAngle = 10.0*kAngTolerance;

  G4ThreeVector nR, nPs, nPe;
  G4ThreeVector norm, sumnorm(0.,0.,0.);

  rho = std::hypot(p.x(),p.y());
  pt  = std::hypot(p.z(),rho-fRtor);

  G4double  distRMax = std::fabs(pt - fRmax);
  if(fRmin) distRMin = std::fabs(pt - fRmin);

  if( rho > delta && pt != 0.0 )
  {
    G4double redFactor= (rho-fRtor)/rho;
    nR = G4ThreeVector( p.x()*redFactor,  // p.x()*(1.-fRtor/rho),
                        p.y()*redFactor,  // p.y()*(1.-fRtor/rho),
                        p.z()          );
    nR *= 1.0/pt;
  }

  if ( fDPhi < twopi ) // && rho ) // old limitation against (0,0,z)
  {
    if ( rho )
    {
      pPhi = std::atan2(p.y(),p.x());

      if(pPhi < fSPhi-delta)            { pPhi += twopi; }
      else if(pPhi > fSPhi+fDPhi+delta) { pPhi -= twopi; }

      distSPhi = std::fabs( pPhi - fSPhi );
      distEPhi = std::fabs(pPhi-fSPhi-fDPhi);
    }
    nPs = G4ThreeVector(std::sin(fSPhi),-std::cos(fSPhi),0);
    nPe = G4ThreeVector(-std::sin(fSPhi+fDPhi),std::cos(fSPhi+fDPhi),0);
  } 
  if( distRMax <= delta )
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  else if( fRmin && (distRMin <= delta) ) // Must not be on both Outer and Inner
  {
    noSurfaces ++;
    sumnorm -= nR;
  }

  //  To be on one of the 'phi' surfaces,
  //  it must be within the 'tube' - with tolerance

  if( (fDPhi < twopi) && (fRmin-delta <= pt) && (pt <= (fRmax+delta)) )
  {
    if (distSPhi <= dAngle)
    {
      noSurfaces ++;
      sumnorm += nPs;
    }
    if (distEPhi <= dAngle) 
    {
      noSurfaces ++;
      sumnorm += nPe;
    }
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
     G4ExceptionDescription ed;
     ed.precision(16);

     EInside  inIt= Inside( p );
     
     if( inIt != kSurface )
     {
        ed << " ERROR>  Surface Normal was called for Torus,"
           << " with point not on surface." << G4endl;
     }
     else
     {
        ed << " ERROR>  Surface Normal has not found a surface, "
           << " despite the point being on the surface. " <<G4endl;
     }

     if( inIt != kInside)
     {
         ed << " Safety (Dist To In)  = " << DistanceToIn(p) << G4endl;
     }
     if( inIt != kOutside)
     {
         ed << " Safety (Dist to Out) = " << DistanceToOut(p) << G4endl;
     }
     ed << " Coordinates of point : " << p << G4endl;
     ed << " Parameters  of solid : " << G4endl << *this << G4endl;

     if( inIt == kSurface )
     {
        G4Exception("G4Torus::SurfaceNormal(p)", "GeomSolids1002",
                    JustWarning, ed,
                    "Failing to find normal, even though point is on surface!");
     }
     else
     {
        static const char* NameInside[3]= { "Inside", "Surface", "Outside" };
        ed << "  The point is " << NameInside[inIt] << " the solid. "<< G4endl;
        G4Exception("G4Torus::SurfaceNormal(p)", "GeomSolids1002",
                    JustWarning, ed, "Point p is not on surface !?" );
     }
#endif
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = sumnorm.unit(); }

  return norm ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Torus::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  ENorm side ;
  G4ThreeVector norm;
  G4double rho,pt,phi;
  G4double distRMin,distRMax,distSPhi,distEPhi,distMin;

  rho = std::hypot(p.x(),p.y());
  pt  = std::hypot(p.z(),rho-fRtor);

#ifdef G4CSGDEBUG
  G4cout << " G4Torus::ApproximateSurfaceNormal called for point " << p
         << G4endl;
#endif
   
  distRMax = std::fabs(pt - fRmax) ;

  if(fRmin)  // First minimum radius
  {
    distRMin = std::fabs(pt - fRmin) ;

    if (distRMin < distRMax)
    {
      distMin = distRMin ;
      side    = kNRMin ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax ;
    }
  }
  else
  {
    distMin = distRMax ;
    side    = kNRMax ;
  }    
  if ( (fDPhi < twopi) && rho )
  {
    phi = std::atan2(p.y(),p.x()) ; // Protected against (0,0,z) (above rho!=0)

    if (phi < 0)  { phi += twopi ; }

    if (fSPhi < 0 )  { distSPhi = std::fabs(phi-(fSPhi+twopi))*rho ; }
    else             { distSPhi = std::fabs(phi-fSPhi)*rho ; }

    distEPhi = std::fabs(phi - fSPhi - fDPhi)*rho ;

    if (distSPhi < distEPhi) // Find new minimum
    {
      if (distSPhi<distMin) side = kNSPhi ;
    }
    else
    {
      if (distEPhi < distMin)  { side = kNEPhi ; }
    }
  }  
  switch (side)
  {
    case kNRMin:      // Inner radius
      norm = G4ThreeVector( -p.x()*(1-fRtor/rho)/pt,
                            -p.y()*(1-fRtor/rho)/pt,
                            -p.z()/pt                 ) ;
      break ;
    case kNRMax:      // Outer radius
      norm = G4ThreeVector( p.x()*(1-fRtor/rho)/pt,
                            p.y()*(1-fRtor/rho)/pt,
                            p.z()/pt                  ) ;
      break;
    case kNSPhi:
      norm = G4ThreeVector(std::sin(fSPhi),-std::cos(fSPhi),0) ;
      break;
    case kNEPhi:
      norm = G4ThreeVector(-std::sin(fSPhi+fDPhi),std::cos(fSPhi+fDPhi),0) ;
      break;
    default:          // Should never reach this case ...
      DumpInfo();
      G4Exception("G4Torus::ApproxSurfaceNormal()",
                  "GeomSolids1002", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break ;
  } 
  return norm ;
}

///////////////////////////////////////////////////////////////////////
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
//        - if valid (phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points

G4double G4Torus::DistanceToIn( const G4ThreeVector& p,
                                const G4ThreeVector& v ) const
{

  G4double snxt=kInfinity, sphi=kInfinity; // snxt = default return value

  G4double  sd[4] ;

  // Precalculated trig for phi intersections - used by r,z intersections to
  //                                            check validity

  G4bool seg;        // true if segmented
  G4double hDPhi;    // half dphi
  G4double cPhi,sinCPhi=0.,cosCPhi=0.;  // central phi

  G4double tolORMin2;  // `generous' radii squared
  G4double tolORMax2;

  G4double Dist,xi,yi,zi,rhoi,it2; // Intersection point variables

  G4double Comp;
  G4double cosSPhi,sinSPhi;       // Trig for phi start intersect
  G4double ePhi,cosEPhi,sinEPhi;  // for phi end intersect

  // Set phi divided flag and precalcs
  //
  if ( fDPhi < twopi )
  {
    seg        = true ;
    hDPhi      = 0.5*fDPhi ;    // half delta phi
    cPhi       = fSPhi + hDPhi ;
    sinCPhi    = std::sin(cPhi) ;
    cosCPhi    = std::cos(cPhi) ;
  }
  else
  {
    seg = false ;
  }

  if (fRmin > fRminTolerance) // Calculate tolerant rmin and rmax
  {
    tolORMin2 = (fRmin - fRminTolerance)*(fRmin - fRminTolerance) ;
  }
  else
  {
    tolORMin2 = 0 ;
  }
  tolORMax2 = (fRmax + fRmaxTolerance)*(fRmax + fRmaxTolerance) ;

  // Intersection with Rmax (possible return) and Rmin (must also check phi)

  snxt = SolveNumericJT(p,v,fRmax,true);

  if (fRmin)  // Possible Rmin intersection
  {
    sd[0] = SolveNumericJT(p,v,fRmin,true);
    if ( sd[0] < snxt )  { snxt = sd[0] ; }
  }

  //
  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to kCarTolerance*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?

  if (seg)
  {
    sinSPhi = std::sin(fSPhi) ; // First phi surface ('S'tarting phi)
    cosSPhi = std::cos(fSPhi) ;
    Comp    = v.x()*sinSPhi - v.y()*cosSPhi ;  // Component in outwards
                                               // normal direction
    if (Comp < 0 )
    {
      Dist = (p.y()*cosSPhi - p.x()*sinSPhi) ;

      if (Dist < halfCarTolerance)
      {
        sphi = Dist/Comp ;
        if (sphi < snxt)
        {
          if ( sphi < 0 )  { sphi = 0 ; }

          xi    = p.x() + sphi*v.x() ;
          yi    = p.y() + sphi*v.y() ;
          zi    = p.z() + sphi*v.z() ;
          rhoi = std::hypot(xi,yi);
          it2 = zi*zi + (rhoi-fRtor)*(rhoi-fRtor);

          if ( it2 >= tolORMin2 && it2 <= tolORMax2 )
          {
            // r intersection is good - check intersecting
            // with correct half-plane
            //
            if ((yi*cosCPhi-xi*sinCPhi)<=0)  { snxt=sphi; }
          }
        }
      }
    }
    ePhi=fSPhi+fDPhi;    // Second phi surface ('E'nding phi)
    sinEPhi=std::sin(ePhi);
    cosEPhi=std::cos(ePhi);
    Comp=-(v.x()*sinEPhi-v.y()*cosEPhi);
        
    if ( Comp < 0 )   // Component in outwards normal dirn
    {
      Dist = -(p.y()*cosEPhi - p.x()*sinEPhi) ;

      if (Dist < halfCarTolerance )
      {
        sphi = Dist/Comp ;

        if (sphi < snxt )
        {
          if (sphi < 0 )  { sphi = 0 ; }
       
          xi    = p.x() + sphi*v.x() ;
          yi    = p.y() + sphi*v.y() ;
          zi    = p.z() + sphi*v.z() ;
          rhoi = std::hypot(xi,yi);
          it2 = zi*zi + (rhoi-fRtor)*(rhoi-fRtor);

          if (it2 >= tolORMin2 && it2 <= tolORMax2)
          {
            // z and r intersections good - check intersecting
            // with correct half-plane
            //
            if ((yi*cosCPhi-xi*sinCPhi)>=0)  { snxt=sphi; }
          }    
        }
      }
    }
  }
  if(snxt < halfCarTolerance)  { snxt = 0.0 ; }

  return snxt ;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Torus::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe=0.0, safe1, safe2 ;
  G4double phiC, cosPhiC, sinPhiC, safePhi, ePhi, cosPsi ;
  G4double rho, pt ;
  
  rho = std::hypot(p.x(),p.y());
  pt  = std::hypot(p.z(),rho-fRtor);
  safe1 = fRmin - pt ;
  safe2 = pt - fRmax ;

  if (safe1 > safe2)  { safe = safe1; }
  else                { safe = safe2; }

  if ( fDPhi < twopi && rho )
  {
    phiC    = fSPhi + fDPhi*0.5 ;
    cosPhiC = std::cos(phiC) ;
    sinPhiC = std::sin(phiC) ;
    cosPsi  = (p.x()*cosPhiC + p.y()*sinPhiC)/rho ;

    if (cosPsi < std::cos(fDPhi*0.5) ) // Psi=angle from central phi to point
    {                                  // Point lies outside phi range
      if ((p.y()*cosPhiC - p.x()*sinPhiC) <= 0 )
      {
        safePhi = std::fabs(p.x()*std::sin(fSPhi) - p.y()*std::cos(fSPhi)) ;
      }
      else
      {
        ePhi    = fSPhi + fDPhi ;
        safePhi = std::fabs(p.x()*std::sin(ePhi) - p.y()*std::cos(ePhi)) ;
      }
      if (safePhi > safe)  { safe = safePhi ; }
    }
  }
  if (safe < 0 )  { safe = 0 ; }
  return safe;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection
//

G4double G4Torus::DistanceToOut( const G4ThreeVector& p,
                                 const G4ThreeVector& v,
                                 const G4bool calcNorm,
                                       G4bool *validNorm,
                                       G4ThreeVector  *n  ) const
{
  ESide    side = kNull, sidephi = kNull ;
  G4double snxt = kInfinity, sphi, sd[4] ;

  // Vars for phi intersection
  //
  G4double sinSPhi, cosSPhi, ePhi, sinEPhi, cosEPhi;
  G4double cPhi, sinCPhi, cosCPhi ;
  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, zi, vphi ;

  // Radial Intersections Defenitions & General Precals

  //////////////////////// new calculation //////////////////////

#if 1

  // This is the version with the calculation of CalcNorm = true 
  // To be done: Check the precision of this calculation.
  // If you want return always validNorm = false, then take the version below
  
  
  G4double rho = std::hypot(p.x(),p.y());
  G4double pt = hypot(p.z(),rho-fRtor);

  G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;

  G4double tolRMax = fRmax - fRmaxTolerance ;
   
  G4double vDotNmax   = pDotV - fRtor*(v.x()*p.x() + v.y()*p.y())/rho ;
  G4double pDotxyNmax = (1 - fRtor/rho) ;

  if( (pt*pt > tolRMax*tolRMax) && (vDotNmax >= 0) )
  {
    // On tolerant boundary & heading outwards (or perpendicular to) outer
    // radial surface -> leaving immediately with *n for really convex part
    // only
      
    if ( calcNorm && (pDotxyNmax >= -2.*fRmaxTolerance) ) 
    {
      *n = G4ThreeVector( p.x()*(1 - fRtor/rho)/pt,
                          p.y()*(1 - fRtor/rho)/pt,
                          p.z()/pt                  ) ;
      *validNorm = true ;
    }
     
    return snxt = 0 ; // Leaving by Rmax immediately
  }
  
  snxt = SolveNumericJT(p,v,fRmax,false);  
  side = kRMax ;

  // rmin

  if ( fRmin )
  {
    G4double tolRMin = fRmin + fRminTolerance ;

    if ( (pt*pt < tolRMin*tolRMin) && (vDotNmax < 0) )
    {
      if (calcNorm)  { *validNorm = false ; } // Concave surface of the torus
      return  snxt = 0 ;                      // Leaving by Rmin immediately
    }
    
    sd[0] = SolveNumericJT(p,v,fRmin,false);
    if ( sd[0] < snxt )
    {
      snxt = sd[0] ;
      side = kRMin ;
    }
  }

#else

  // this is the "conservative" version which return always validnorm = false
  // NOTE: using this version the unit test testG4Torus will break

  snxt = SolveNumericJT(p,v,fRmax,false);  
  side = kRMax ;

  if ( fRmin )
  {
    sd[0] = SolveNumericJT(p,v,fRmin,false);
    if ( sd[0] < snxt )
    {
      snxt = sd[0] ;
      side = kRMin ;
    }
  }

  if ( calcNorm && (snxt == 0.0) )
  {
    *validNorm = false ;    // Leaving solid, but possible re-intersection
    return snxt  ;
  }

#endif
  
  if (fDPhi < twopi)  // Phi Intersections
  {
    sinSPhi = std::sin(fSPhi) ;
    cosSPhi = std::cos(fSPhi) ;
    ePhi    = fSPhi + fDPhi ;
    sinEPhi = std::sin(ePhi) ;
    cosEPhi = std::cos(ePhi) ;
    cPhi    = fSPhi + fDPhi*0.5 ;
    sinCPhi = std::sin(cPhi) ;
    cosCPhi = std::cos(cPhi) ;
    
    // angle calculation with correction 
    // of difference in domain of atan2 and Sphi
    //
    vphi = std::atan2(v.y(),v.x()) ;
     
    if ( vphi < fSPhi - halfAngTolerance  )    { vphi += twopi; }
    else if ( vphi > ePhi + halfAngTolerance ) { vphi -= twopi; }

    if ( p.x() || p.y() ) // Check if on z axis (rho not needed later)
    {
      pDistS = p.x()*sinSPhi - p.y()*cosSPhi ; // pDist -ve when inside
      pDistE = -p.x()*sinEPhi + p.y()*cosEPhi ;

      // Comp -ve when in direction of outwards normal
      //
      compS   = -sinSPhi*v.x() + cosSPhi*v.y() ;
      compE   = sinEPhi*v.x() - cosEPhi*v.y() ;
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
              sidephi = kSPhi;
              if ( ((fSPhi-halfAngTolerance)<=vphi)
                && ((ePhi+halfAngTolerance)>=vphi) )
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
          if ( (sphi2 > -kCarTolerance) && (sphi2 < sphi) )
          {
            xi = p.x() + sphi2*v.x() ;
            yi = p.y() + sphi2*v.y() ;
              
            if ( (std::fabs(xi)<=kCarTolerance)
              && (std::fabs(yi)<=kCarTolerance) )
            {
              // Leaving via ending phi
              //
              if( !( (fSPhi-halfAngTolerance <= vphi)
                  && (ePhi+halfAngTolerance >= vphi) ) )
              {
                sidephi = kEPhi ;
                sphi = sphi2;
              }
            } 
            else    // Check intersecting with correct half-plane 
            {
              if ( (yi*cosCPhi-xi*sinCPhi) >= 0)
              {
                // Leaving via ending phi
                //
                sidephi = kEPhi ;
                sphi = sphi2;
               
              }
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

      vphi = std::atan2(v.y(),v.x());
 
      if ( ( fSPhi-halfAngTolerance <= vphi ) && 
           ( vphi <= ( ePhi+halfAngTolerance ) ) )
      {
        sphi = kInfinity;
      }
      else
      {
        sidephi = kSPhi ; // arbitrary 
        sphi=0;
      }
    }

    // Order intersections

    if (sphi<snxt)
    {
      snxt=sphi;
      side=sidephi;
    }     
  }

  G4double rhoi,it,iDotxyNmax ;
  // Note: by numerical computation we know where the ray hits the torus
  // So I propose to return the side where the ray hits

  if (calcNorm)
  {
    switch(side)
    {
      case kRMax:                     // n is unit vector 
        xi    = p.x() + snxt*v.x() ;
        yi    = p.y() + snxt*v.y() ;
        zi    = p.z() + snxt*v.z() ;
        rhoi = std::hypot(xi,yi);
        it = hypot(zi,rhoi-fRtor);

        iDotxyNmax = (1-fRtor/rhoi) ;
        if(iDotxyNmax >= -2.*fRmaxTolerance) // really convex part of Rmax
        {                       
          *n = G4ThreeVector( xi*(1-fRtor/rhoi)/it,
                              yi*(1-fRtor/rhoi)/it,
                              zi/it                 ) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ; // concave-convex part of Rmax
        }
        break ;

      case kRMin:
        *validNorm = false ;  // Rmin is concave or concave-convex
        break;

      case kSPhi:
        if (fDPhi <= pi )
        {
          *n=G4ThreeVector(std::sin(fSPhi),-std::cos(fSPhi),0);
          *validNorm=true;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kEPhi:
        if (fDPhi <= pi)
        {
          *n=G4ThreeVector(-std::sin(fSPhi+fDPhi),std::cos(fSPhi+fDPhi),0);
          *validNorm=true;
        }
        else
        {
          *validNorm = false ;
        }
        break;

      default:

        // It seems we go here from time to time ...

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
                << "snxt = " << snxt/mm << " mm" << G4endl;
        message.precision(oldprc);
        G4Exception("G4Torus::DistanceToOut(p,v,..)",
                    "GeomSolids1002",JustWarning, message);
        break;
    }
  }
  if ( snxt<halfCarTolerance )  { snxt=0 ; }

  return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Torus::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0,safeR1,safeR2;
  G4double rho,pt ;
  G4double safePhi,phiC,cosPhiC,sinPhiC,ePhi;
  
  rho = std::hypot(p.x(),p.y());
  pt  = std::hypot(p.z(),rho-fRtor);
  
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
     G4cout.precision(oldprc);
     G4Exception("G4Torus::DistanceToOut(p)", "GeomSolids1002",
                 JustWarning, "Point p is outside !?" );
  }
#endif

  if (fRmin)
  {
    safeR1 = pt - fRmin ;
    safeR2 = fRmax - pt ;

    if (safeR1 < safeR2)  { safe = safeR1 ; }
    else                  { safe = safeR2 ; }
  }
  else
  {
    safe = fRmax - pt ;
  }  

  // Check if phi divided, Calc distances closest phi plane
  //
  if (fDPhi<twopi) // Above/below central phi of Torus?
  {
    phiC    = fSPhi + fDPhi*0.5 ;
    cosPhiC = std::cos(phiC) ;
    sinPhiC = std::sin(phiC) ;

    if ((p.y()*cosPhiC-p.x()*sinPhiC)<=0)
    {
      safePhi = -(p.x()*std::sin(fSPhi) - p.y()*std::cos(fSPhi)) ;
    }
    else
    {
      ePhi    = fSPhi + fDPhi ;
      safePhi = (p.x()*std::sin(ePhi) - p.y()*std::cos(ePhi)) ;
    }
    if (safePhi < safe)  { safe = safePhi ; }
  }
  if (safe < 0)  { safe = 0 ; }
  return safe ;  
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

G4GeometryType G4Torus::GetEntityType() const
{
  return G4String("G4Torus");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Torus::Clone() const
{
  return new G4Torus(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Torus::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Torus\n"
     << " Parameters: \n"
     << "    inner radius: " << fRmin/mm << " mm \n"
     << "    outer radius: " << fRmax/mm << " mm \n"
     << "    swept radius: " << fRtor/mm << " mm \n"
     << "    starting phi: " << fSPhi/degree << " degrees \n"
     << "    delta phi   : " << fDPhi/degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Torus::GetPointOnSurface() const
{
  G4double cosu, sinu,cosv, sinv, aOut, aIn, aSide, chose, phi, theta, rRand;
   
  phi   = G4RandFlat::shoot(fSPhi,fSPhi+fDPhi);
  theta = G4RandFlat::shoot(0.,twopi);
  
  cosu   = std::cos(phi);    sinu = std::sin(phi);
  cosv   = std::cos(theta);  sinv = std::sin(theta); 

  // compute the areas

  aOut   = (fDPhi)*twopi*fRtor*fRmax;
  aIn    = (fDPhi)*twopi*fRtor*fRmin;
  aSide  = pi*(fRmax*fRmax-fRmin*fRmin);
  
  if ((fSPhi == 0) && (fDPhi == twopi)){ aSide = 0; }
  chose = G4RandFlat::shoot(0.,aOut + aIn + 2.*aSide);

  if(chose < aOut)
  {
    return G4ThreeVector ((fRtor+fRmax*cosv)*cosu,
                          (fRtor+fRmax*cosv)*sinu, fRmax*sinv);
  }
  else if( (chose >= aOut) && (chose < aOut + aIn) )
  {
    return G4ThreeVector ((fRtor+fRmin*cosv)*cosu,
                          (fRtor+fRmin*cosv)*sinu, fRmin*sinv);
  }
  else if( (chose >= aOut + aIn) && (chose < aOut + aIn + aSide) )
  {
    rRand = GetRadiusInRing(fRmin,fRmax);
    return G4ThreeVector ((fRtor+rRand*cosv)*std::cos(fSPhi),
                          (fRtor+rRand*cosv)*std::sin(fSPhi), rRand*sinv);
  }
  else
  {   
    rRand = GetRadiusInRing(fRmin,fRmax);
    return G4ThreeVector ((fRtor+rRand*cosv)*std::cos(fSPhi+fDPhi),
                          (fRtor+rRand*cosv)*std::sin(fSPhi+fDPhi), 
                          rRand*sinv);
   }
}

///////////////////////////////////////////////////////////////////////
//
// Visualisation Functions

void G4Torus::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Torus::CreatePolyhedron () const 
{
  return new G4PolyhedronTorus (fRmin, fRmax, fRtor, fSPhi, fDPhi);
}

#endif // !defined(G4GEOM_USE_TORUS) || !defined(G4GEOM_USE_SYS_USOLIDS)
