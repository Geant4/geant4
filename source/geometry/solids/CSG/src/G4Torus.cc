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
// $Id: G4Torus.cc,v 1.71 2010-10-19 15:42:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Torus
//
// Implementation
//
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

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4JTPolynomialSolver.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBStube.hh"
#include "G4NURBScylinder.hh"
#include "G4NURBStubesector.hh"

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
  fpPolyhedron = 0;

  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  
  if ( pRtor >= pRmax+1.e3*kCarTolerance )  // Check swept radius, as in G4Cons
  {
    fRtor = pRtor ;
  }
  else
  {
    G4cerr << "ERROR - G4Torus()::SetAllParameters(): " << GetName() << G4endl
           << "        Invalid swept radius !" << G4endl
           << "pRtor = " << pRtor << ", pRmax = " << pRmax << G4endl;
    G4Exception("G4Torus::SetAllParameters()",
                "InvalidSetup", FatalException, "Invalid swept radius.");
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
    G4cerr << "ERROR - G4Torus()::SetAllParameters(): " << GetName() << G4endl
           << "        Invalid values for radii !" << G4endl
           << "        pRmin = " << pRmin << ", pRmax = " << pRmax << G4endl;
    G4Exception("G4Torus::SetAllParameters()",
                "InvalidSetup", FatalException, "Invalid radii.");
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
      G4cerr << "ERROR - G4Torus::SetAllParameters(): " << GetName() << G4endl
             << "        Negative Z delta-Phi ! - "
             << pDPhi << G4endl;
      G4Exception("G4Torus::SetAllParameters()",
                  "InvalidSetup", FatalException, "Invalid dphi.");
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
    kRadTolerance(0.), kAngTolerance(0.)
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
    kRadTolerance(rhs.kRadTolerance), kAngTolerance(rhs.kAngTolerance)
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
  G4double c[5], sr[4], si[4] ;

  G4double Rtor2 = fRtor*fRtor, r2 = r*r  ;

  G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
  G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;

  c[0] = 1.0 ;
  c[1] = 4*pDotV ;
  c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - r2 + 2*Rtor2*v.z()*v.z()) ;
  c[3] = 4*(pDotV*(pRad2 - Rtor2 - r2) + 2*Rtor2*p.z()*v.z()) ;
  c[4] = pRad2*pRad2 - 2*pRad2*(Rtor2+r2) 
       + 4*Rtor2*p.z()*p.z() + (Rtor2-r2)*(Rtor2-r2) ;
  
  G4JTPolynomialSolver  torusEq;

  num = torusEq.FindRoots( c, 4, sr, si );
  
  for ( i = 0; i < num; i++ ) 
  {
    if( si[i] == 0. )  { roots.push_back(sr[i]) ; }  // store real roots
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
  static const G4double halfCarTolerance = 0.5*kCarTolerance;
  static const G4double halfAngTolerance = 0.5*kAngTolerance;

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

          scal = v* G4ThreeVector( p.x()*(1-fRtor/std::sqrt(p.x()*p.x()
                                          + p.y()*p.y())),
                                   p.y()*(1-fRtor/std::sqrt(p.x()*p.x()
                                          + p.y()*p.y())),
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
          scal = v* G4ThreeVector( p.x()*(1-fRtor/std::sqrt(p.x()*p.x()
                                          + p.y()*p.y())),
                                   p.y()*(1-fRtor/std::sqrt(p.x()*p.x()
                                          + p.y()*p.y())),
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
// Calculate extent under transform and specified limit

G4bool G4Torus::CalculateExtent( const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                       G4double& pMin, G4double& pMax) const
{
  if ((!pTransform.IsRotated()) && (fDPhi==twopi) && (fRmin==0))
  {
    // Special case handling for unrotated solid torus
    // Compute x/y/z mins and maxs for bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case
      
    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;

    G4double RTorus,delta,diff1,diff2,maxDiff,newMin,newMax;
    G4double xoff1,xoff2,yoff1,yoff2;

    xoffset = pTransform.NetTranslation().x();
    xMin    = xoffset - fRmax - fRtor ;
    xMax    = xoffset + fRmax + fRtor ;

    if (pVoxelLimit.IsXLimited())
    {
      if ( (xMin > pVoxelLimit.GetMaxXExtent()+kCarTolerance)
        || (xMax < pVoxelLimit.GetMinXExtent()-kCarTolerance) )
        return false ;
      else
      {
        if (xMin < pVoxelLimit.GetMinXExtent())
        {
          xMin = pVoxelLimit.GetMinXExtent() ;
        }
        if (xMax > pVoxelLimit.GetMaxXExtent())
        {
          xMax = pVoxelLimit.GetMaxXExtent() ;
        }
      }
    }
    yoffset = pTransform.NetTranslation().y();
    yMin    = yoffset - fRmax - fRtor ;
    yMax    = yoffset + fRmax + fRtor ;

    if (pVoxelLimit.IsYLimited())
    {
      if ( (yMin > pVoxelLimit.GetMaxYExtent()+kCarTolerance)
        || (yMax < pVoxelLimit.GetMinYExtent()-kCarTolerance) )
      {
        return false ;
      }
      else
      {
        if (yMin < pVoxelLimit.GetMinYExtent() )
        { 
          yMin = pVoxelLimit.GetMinYExtent() ;
        }
        if (yMax > pVoxelLimit.GetMaxYExtent() )
        {
          yMax = pVoxelLimit.GetMaxYExtent() ;
        }
      }
    }
    zoffset = pTransform.NetTranslation().z() ;
    zMin    = zoffset - fRmax ;
    zMax    = zoffset + fRmax ;

    if (pVoxelLimit.IsZLimited())
    {
      if ( (zMin > pVoxelLimit.GetMaxZExtent()+kCarTolerance)
        || (zMax < pVoxelLimit.GetMinZExtent()-kCarTolerance) )
      {
        return false ;
      }
      else
      {
        if (zMin < pVoxelLimit.GetMinZExtent() )
        {
          zMin = pVoxelLimit.GetMinZExtent() ;
        }
        if (zMax > pVoxelLimit.GetMaxZExtent() )
        {
          zMax = pVoxelLimit.GetMaxZExtent() ;
        }
      }
    }

    // Known to cut cylinder
    
    switch (pAxis)
    {
      case kXAxis:
        yoff1=yoffset-yMin;
        yoff2=yMax-yoffset;
        if ( yoff1 >= 0 && yoff2 >= 0 )
        {
          // Y limits cross max/min x => no change
          //
          pMin = xMin ;
          pMax = xMax ;
        }
        else
        {
          // Y limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs
          //
       
          RTorus=fRmax+fRtor;
          delta   = RTorus*RTorus - yoff1*yoff1;
          diff1   = (delta>0.) ? std::sqrt(delta) : 0.;
          delta   = RTorus*RTorus - yoff2*yoff2;
          diff2   = (delta>0.) ? std::sqrt(delta) : 0.;
          maxDiff = (diff1 > diff2) ? diff1:diff2 ;
          newMin  = xoffset - maxDiff ;
          newMax  = xoffset + maxDiff ;
          pMin    = (newMin < xMin) ? xMin : newMin ;
          pMax    = (newMax > xMax) ? xMax : newMax ;
        }
        break;

      case kYAxis:
        xoff1 = xoffset - xMin ;
        xoff2 = xMax - xoffset ;
        if (xoff1 >= 0 && xoff2 >= 0 )
        {
          // X limits cross max/min y => no change
          //
          pMin = yMin ;
          pMax = yMax ;
        } 
        else
        {
          // X limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs
          //
          RTorus=fRmax+fRtor;
          delta   = RTorus*RTorus - xoff1*xoff1;
          diff1   = (delta>0.) ? std::sqrt(delta) : 0.;
          delta   = RTorus*RTorus - xoff2*xoff2;
          diff2   = (delta>0.) ? std::sqrt(delta) : 0.;
          maxDiff = (diff1 > diff2) ? diff1 : diff2 ;
          newMin  = yoffset - maxDiff ;
          newMax  = yoffset + maxDiff ;
          pMin    = (newMin < yMin) ? yMin : newMin ;
          pMax    = (newMax > yMax) ? yMax : newMax ;
        }
        break;

      case kZAxis:
        pMin=zMin;
        pMax=zMax;
        break;

      default:
        break;
    }
    pMin -= kCarTolerance ;
    pMax += kCarTolerance ;

    return true;
  }
  else
  {
    G4int i, noEntries, noBetweenSections4 ;
    G4bool existsAfterClip = false ;

    // Calculate rotated vertex coordinates

    G4ThreeVectorList *vertices ;
    G4int noPolygonVertices ;  // will be 4 
    vertices = CreateRotatedVertices(pTransform,noPolygonVertices) ;

    pMin = +kInfinity ;
    pMax = -kInfinity ;

    noEntries          = vertices->size() ;
    noBetweenSections4 = noEntries - noPolygonVertices ;
      
    for (i=0;i<noEntries;i+=noPolygonVertices)
    {
      ClipCrossSection(vertices,i,pVoxelLimit,pAxis,pMin,pMax);
    }    
    for (i=0;i<noBetweenSections4;i+=noPolygonVertices)
    {
      ClipBetweenSections(vertices,i,pVoxelLimit,pAxis,pMin,pMax);
    }
    if (pMin!=kInfinity||pMax!=-kInfinity)
    {
      existsAfterClip = true ; // Add 2*tolerance to avoid precision troubles
      pMin           -= kCarTolerance ;
      pMax           += kCarTolerance ;
    }
    else
    {
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.

      G4ThreeVector clipCentre(
          (pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
          (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
          (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5  ) ;
        
      if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside )
      {
        existsAfterClip = true ;
        pMin            = pVoxelLimit.GetMinExtent(pAxis) ;
        pMax            = pVoxelLimit.GetMaxExtent(pAxis) ;
      }
    }
    delete vertices;
    return existsAfterClip;
  }
}

//////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Torus::Inside( const G4ThreeVector& p ) const
{
  G4double r2, pt2, pPhi, tolRMin, tolRMax ;

  EInside in = kOutside ;
  static const G4double halfAngTolerance = 0.5*kAngTolerance;

                                               // General precals
  r2  = p.x()*p.x() + p.y()*p.y() ;
  pt2 = r2 + p.z()*p.z() + fRtor*fRtor - 2*fRtor*std::sqrt(r2) ;

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
  G4double rho2, rho, pt2, pt, pPhi;
  G4double distRMin = kInfinity;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;

  static const G4double delta = 0.5*kCarTolerance;
  static const G4double dAngle = 0.5*kAngTolerance;

  G4ThreeVector nR, nPs, nPe;
  G4ThreeVector norm, sumnorm(0.,0.,0.);

  rho2 = p.x()*p.x() + p.y()*p.y();
  rho = std::sqrt(rho2);
  pt2 = std::fabs(rho2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*rho);
  pt = std::sqrt(pt2) ;

  G4double  distRMax = std::fabs(pt - fRmax);
  if(fRmin) distRMin = std::fabs(pt - fRmin);

  if( rho > delta )
  {
    nR = G4ThreeVector( p.x()*(1-fRtor/rho)/pt,
                        p.y()*(1-fRtor/rho)/pt,
                        p.z()/pt                 );
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
  if( fRmin && distRMin <= delta )
  {
    noSurfaces ++;
    sumnorm -= nR;
  }
  if( fDPhi < twopi )   
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
    G4Exception("G4Torus::SurfaceNormal(p)", "Notification", JustWarning, 
                "Point p is not on surface !?" );
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
  G4double rho2,rho,pt2,pt,phi;
  G4double distRMin,distRMax,distSPhi,distEPhi,distMin;

  rho2 = p.x()*p.x() + p.y()*p.y();
  rho = std::sqrt(rho2) ;
  pt2 = std::fabs(rho2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*rho) ;
  pt = std::sqrt(pt2) ;

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
                  "Notification", JustWarning,
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

  G4double  s[4] ;

  // Precalculated trig for phi intersections - used by r,z intersections to
  //                                            check validity

  G4bool seg;        // true if segmented
  G4double hDPhi;    // half dphi
  G4double cPhi,sinCPhi=0.,cosCPhi=0.;  // central phi

  G4double tolORMin2;  // `generous' radii squared
  G4double tolORMax2;

  static const G4double halfCarTolerance = 0.5*kCarTolerance;
 
  G4double Dist,xi,yi,zi,rhoi2,it2; // Intersection point variables

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

  G4double Rtor2 = fRtor*fRtor ;

  snxt = SolveNumericJT(p,v,fRmax,true);
  if (fRmin)  // Possible Rmin intersection
  {
    s[0] = SolveNumericJT(p,v,fRmin,true);
    if ( s[0] < snxt )  { snxt = s[0] ; }
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
    sinSPhi = std::sin(fSPhi) ; // First phi surface (`S'tarting phi)
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
          rhoi2 = xi*xi + yi*yi ;
          it2   = std::fabs(rhoi2 + zi*zi + Rtor2 - 2*fRtor*std::sqrt(rhoi2)) ;

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
    ePhi=fSPhi+fDPhi;    // Second phi surface (`E'nding phi)
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
          rhoi2 = xi*xi + yi*yi ;
          it2   = std::fabs(rhoi2 + zi*zi + Rtor2 - 2*fRtor*std::sqrt(rhoi2)) ;

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
  G4double rho2, rho, pt2, pt ;
    
  rho2 = p.x()*p.x() + p.y()*p.y() ;
  rho  = std::sqrt(rho2) ;
  pt2  = std::fabs(rho2 + p.z()*p.z() + fRtor*fRtor - 2*fRtor*rho) ;
  pt   = std::sqrt(pt2) ;

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
  G4double snxt = kInfinity, sphi, s[4] ;

  static const G4double halfCarTolerance = 0.5*kCarTolerance;
  static const G4double halfAngTolerance = 0.5*kAngTolerance;

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
  
  G4double Rtor2 = fRtor*fRtor ;
  G4double rho2  = p.x()*p.x()+p.y()*p.y();
  G4double rho   = std::sqrt(rho2) ;


  G4double pt2   = std::fabs(rho2 + p.z()*p.z() + Rtor2 - 2*fRtor*rho) ;
  G4double pt    = std::sqrt(pt2) ;

  G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;

  G4double tolRMax = fRmax - fRmaxTolerance ;
   
  G4double vDotNmax   = pDotV - fRtor*(v.x()*p.x() + v.y()*p.y())/rho ;
  G4double pDotxyNmax = (1 - fRtor/rho) ;

  if( (pt2 > tolRMax*tolRMax) && (vDotNmax >= 0) )
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

    if ( (pt2 < tolRMin*tolRMin) && (vDotNmax < 0) )
    {
      if (calcNorm)  { *validNorm = false ; } // Concave surface of the torus
      return  snxt = 0 ;                      // Leaving by Rmin immediately
    }
    
    s[0] = SolveNumericJT(p,v,fRmin,false);
    if ( s[0] < snxt )
    {
      snxt = s[0] ;
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
    s[0] = SolveNumericJT(p,v,fRmin,false);
    if ( s[0] < snxt )
    {
      snxt = s[0] ;
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

  G4double rhoi2,rhoi,it2,it,iDotxyNmax ;
  // Note: by numerical computation we know where the ray hits the torus
  // So I propose to return the side where the ray hits

  if (calcNorm)
  {
    switch(side)
    {
      case kRMax:                     // n is unit vector 
        xi    = p.x() + snxt*v.x() ;
        yi    =p.y() + snxt*v.y() ;
        zi    = p.z() + snxt*v.z() ;
        rhoi2 = xi*xi + yi*yi ;
        rhoi  = std::sqrt(rhoi2) ;
        it2   = std::fabs(rhoi2 + zi*zi + fRtor*fRtor - 2*fRtor*rhoi) ;
        it    = std::sqrt(it2) ;
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

        G4int oldprc = G4cout.precision(16);
        G4cout << G4endl;
        DumpInfo();
        G4cout << "Position:"  << G4endl << G4endl;
        G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
        G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
        G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
        G4cout << "Direction:" << G4endl << G4endl;
        G4cout << "v.x() = "   << v.x() << G4endl;
        G4cout << "v.y() = "   << v.y() << G4endl;
        G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
        G4cout << "Proposed distance :" << G4endl << G4endl;
        G4cout << "snxt = " << snxt/mm << " mm" << G4endl << G4endl;
        G4cout.precision(oldprc);
        G4Exception("G4Torus::DistanceToOut(p,v,..)",
                    "Notification",JustWarning,
                    "Undefined side for valid surface normal to solid.");
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
  G4double rho2,rho,pt2,pt ;
  G4double safePhi,phiC,cosPhiC,sinPhiC,ePhi;
  rho2 = p.x()*p.x() + p.y()*p.y() ;
  rho  = std::sqrt(rho2) ;
  pt2  = std::fabs(rho2 + p.z()*p.z() + fRtor*fRtor - 2*fRtor*rho) ;
  pt   = std::sqrt(pt2) ;

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
     G4Exception("G4Torus::DistanceToOut(p)", "Notification",
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

/////////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fRtor cross section
//          [4-7] +fRtor cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

G4ThreeVectorList*
G4Torus::CreateRotatedVertices( const G4AffineTransform& pTransform,
                                      G4int& noPolygonVertices ) const
{
  G4ThreeVectorList *vertices;
  G4ThreeVector vertex0,vertex1,vertex2,vertex3;
  G4double meshAngle,meshRMax,crossAngle,cosCrossAngle,sinCrossAngle,sAngle;
  G4double rMaxX,rMaxY,rMinX,rMinY;
  G4int crossSection,noCrossSections;

  // Compute no of cross-sections necessary to mesh tube
  //
  noCrossSections = G4int (fDPhi/kMeshAngleDefault) + 1 ;

  if (noCrossSections < kMinMeshSections)
  {
    noCrossSections = kMinMeshSections ;
  }
  else if (noCrossSections>kMaxMeshSections)
  {
    noCrossSections=kMaxMeshSections;
  }
  meshAngle = fDPhi/(noCrossSections - 1) ;
  meshRMax  = (fRtor + fRmax)/std::cos(meshAngle*0.5) ;

  // If complete in phi, set start angle such that mesh will be at fRmax
  // on the x axis. Will give better extent calculations when not rotated

  if ( (fDPhi == twopi) && (fSPhi == 0) )
  {
    sAngle = -meshAngle*0.5 ;
  }
  else
  {
    sAngle = fSPhi ;
  }
  vertices = new G4ThreeVectorList();
  
  if (vertices)
  {
    vertices->reserve(noCrossSections*4) ;
    for (crossSection=0;crossSection<noCrossSections;crossSection++)
    {
      // Compute coordinates of cross section at section crossSection

      crossAngle=sAngle+crossSection*meshAngle;
      cosCrossAngle=std::cos(crossAngle);
      sinCrossAngle=std::sin(crossAngle);

      rMaxX=meshRMax*cosCrossAngle;
      rMaxY=meshRMax*sinCrossAngle;
      rMinX=(fRtor-fRmax)*cosCrossAngle;
      rMinY=(fRtor-fRmax)*sinCrossAngle;
      vertex0=G4ThreeVector(rMinX,rMinY,-fRmax);
      vertex1=G4ThreeVector(rMaxX,rMaxY,-fRmax);
      vertex2=G4ThreeVector(rMaxX,rMaxY,+fRmax);
      vertex3=G4ThreeVector(rMinX,rMinY,+fRmax);

      vertices->push_back(pTransform.TransformPoint(vertex0));
      vertices->push_back(pTransform.TransformPoint(vertex1));
      vertices->push_back(pTransform.TransformPoint(vertex2));
      vertices->push_back(pTransform.TransformPoint(vertex3));
    }
    noPolygonVertices = 4 ;
  }
  else
  {
    DumpInfo();
    G4Exception("G4Torus::CreateRotatedVertices()",
                "FatalError", FatalException,
                "Error in allocation of vertices. Out of memory !");
  }
  return vertices;
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

  return os;
}

////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Torus::GetPointOnSurface() const
{
  G4double cosu, sinu,cosv, sinv, aOut, aIn, aSide, chose, phi, theta, rRand;
   
  phi   = RandFlat::shoot(fSPhi,fSPhi+fDPhi);
  theta = RandFlat::shoot(0.,twopi);
  
  cosu   = std::cos(phi);    sinu = std::sin(phi);
  cosv   = std::cos(theta);  sinv = std::sin(theta); 

  // compute the areas

  aOut   = (fDPhi)*twopi*fRtor*fRmax;
  aIn    = (fDPhi)*twopi*fRtor*fRmin;
  aSide  = pi*(fRmax*fRmax-fRmin*fRmin);
  
  if ((fSPhi == 0) && (fDPhi == twopi)){ aSide = 0; }
  chose = RandFlat::shoot(0.,aOut + aIn + 2.*aSide);

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
    rRand = RandFlat::shoot(fRmin,fRmax);
    return G4ThreeVector ((fRtor+rRand*cosv)*std::cos(fSPhi),
                          (fRtor+rRand*cosv)*std::sin(fSPhi), rRand*sinv);
  }
  else
  {   
    rRand = RandFlat::shoot(fRmin,fRmax);
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

G4NURBS* G4Torus::CreateNURBS () const 
{
  G4NURBS* pNURBS;
  if (fRmin != 0) 
  {
    if (fDPhi >= twopi) 
    {
      pNURBS = new G4NURBStube(fRmin, fRmax, fRtor);
    }
    else 
    {
      pNURBS = new G4NURBStubesector(fRmin, fRmax, fRtor, fSPhi, fSPhi + fDPhi);
    }
  }
  else 
  {
    if (fDPhi >= twopi) 
    {
      pNURBS = new G4NURBScylinder (fRmax, fRtor);
    }
    else 
    {
      const G4double epsilon = 1.e-4; // Cylinder sector not yet available!
      pNURBS = new G4NURBStubesector (epsilon, fRmax, fRtor,
                                      fSPhi, fSPhi + fDPhi);
    }
  }
  return pNURBS;
}
