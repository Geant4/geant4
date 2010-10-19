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
// $Id: G4Para.cc,v 1.43 2010-10-19 15:42:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Para
//
// Implementation for G4Para class
//
// History:
//
// 23.10.05 V.Grichine: bug fixed in DistanceToOut(p,v,...) for the v.x()<0 case 
// 28.04.05 V.Grichine: new SurfaceNormal according to J. Apostolakis proposal 
// 30.11.04 V.Grichine: modifications in SurfaceNormal for edges/vertices and
//                      in constructor with vertices
// 14.02.02 V.Grichine: bug fixed in Inside according to proposal of D.Wright
// 18.11.99 V.Grichine: kUndef was added to ESide
// 31.10.96 V.Grichine: Modifications according G4Box/Tubs before to commit
// 21.03.95 P.Kent: Modified for `tolerant' geom
//
////////////////////////////////////////////////////////////////////////////

#include "G4Para.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"

using namespace CLHEP;

// Private enum: Not for external use 
    
enum ESide {kUndef,kPX,kMX,kPY,kMY,kPZ,kMZ};

// used internally for normal routine

enum ENSide {kNZ,kNX,kNY};

/////////////////////////////////////////////////////////////////////
//
// Constructor - check and set half-widths

void G4Para::SetAllParameters( G4double pDx, G4double pDy, G4double pDz, 
                               G4double pAlpha, G4double pTheta, G4double pPhi )
{
  if ( pDx > 0 && pDy > 0 && pDz > 0 )
  {
    fDx         = pDx;
    fDy         = pDy;
    fDz         = pDz;
    fTalpha     = std::tan(pAlpha);
    fTthetaCphi = std::tan(pTheta)*std::cos(pPhi);
    fTthetaSphi = std::tan(pTheta)*std::sin(pPhi);
  }
  else
  {
    G4cerr << "ERROR - G4Para()::SetAllParameters(): " << GetName() << G4endl
           << "        Invalid dimensions ! - "
           << pDx << ", " << pDy << ", " << pDz << G4endl;
    G4Exception("G4Para::SetAllParameters()", "InvalidSetup",
                FatalException, "Invalid Length Parameters.");
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////////////
//

G4Para::G4Para(const G4String& pName,
                     G4double pDx, G4double pDy, G4double pDz,
                     G4double pAlpha, G4double pTheta, G4double pPhi)
  : G4CSGSolid(pName)
{
  if ((pDx<=0) || (pDy<=0) || (pDz<=0))
  {
    G4cerr << "ERROR - G4Para()::G4Para(): " << GetName() << G4endl
           << "        Invalid dimensions ! - "
           << pDx << ", " << pDy << ", " << pDz << G4endl;
    G4Exception("G4Para::G4Para()", "InvalidSetup",
                FatalException, "Invalid Length Parameters.");
  }
  SetAllParameters( pDx, pDy, pDz, pAlpha, pTheta, pPhi);
}

////////////////////////////////////////////////////////////////////////
//
// Constructor - Design of trapezoid based on 8 G4ThreeVector parameters, 
// which are its vertices. Checking of planarity with preparation of 
// fPlanes[] and than calculation of other members

G4Para::G4Para( const G4String& pName,
                const G4ThreeVector pt[8] )
  : G4CSGSolid(pName)
{
  if (!( pt[0].z()<0 && pt[0].z()==pt[1].z() && pt[0].z()==pt[2].z() &&
         pt[0].z()==pt[3].z() && pt[4].z()>0 && pt[4].z()==pt[5].z() &&
         pt[4].z()==pt[6].z() && pt[4].z()==pt[7].z()           &&
        (pt[0].z()+pt[4].z())==0                                &&
         pt[0].y()==pt[1].y() && pt[2].y()==pt[3].y()           &&
         pt[4].y()==pt[5].y() && pt[6].y()==pt[7].y()           &&
       ( pt[0].y() + pt[2].y() + pt[4].y() + pt[6].y() ) == 0   && 
       ( pt[0].x() + pt[1].x() + pt[4].x() + pt[5].x() ) == 0) )
  {
    G4cerr << "ERROR - G4Para()::G4Para(): " << GetName() << G4endl
           << "        Invalid dimensions !" << G4endl;
    G4Exception("G4Para::G4Para()", "InvalidSetup",
                FatalException, "Invalid vertice coordinates.");
  }    
  fDx = ((pt[3]).x()-(pt[2]).x())*0.5;
  fDy = ((pt[2]).y()-(pt[1]).y())*0.5;
  fDz = (pt[7]).z();
  fTalpha = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/fDy ;
  fTthetaCphi = ((pt[4]).x()+fDy*fTalpha+fDx)/fDz ;
  fTthetaSphi = ((pt[4]).y()+fDy)/fDz ;
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Para::G4Para( __void__& a )
  : G4CSGSolid(a), fDx(0.), fDy(0.), fDz(0.),
    fTalpha(0.), fTthetaCphi(0.), fTthetaSphi(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//

G4Para::~G4Para()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Para::G4Para(const G4Para& rhs)
  : G4CSGSolid(rhs), fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz),
    fTalpha(rhs.fTalpha), fTthetaCphi(rhs.fTthetaCphi),
    fTthetaSphi(rhs.fTthetaSphi)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Para& G4Para::operator = (const G4Para& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fDx = rhs.fDx; fDy = rhs.fDy; fDz = rhs.fDz;
   fTalpha = rhs.fTalpha; fTthetaCphi = rhs.fTthetaCphi;
   fTthetaSphi = rhs.fTthetaSphi;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Para::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}


//////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Para::CalculateExtent( const EAxis pAxis,
                                const G4VoxelLimits& pVoxelLimit,
                                const G4AffineTransform& pTransform,
                                     G4double& pMin, G4double& pMax ) const
{
  G4bool flag;

  if (!pTransform.IsRotated())
  {  
    // Special case handling for unrotated trapezoids
    // Compute z/x/y/ mins and maxs respecting limits, with early returns
    // if outside limits. Then switch() on pAxis

    G4int i ; 
    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;
    G4double temp[8] ;       // some points for intersection with zMin/zMax
    
    xoffset=pTransform.NetTranslation().x();      
    yoffset=pTransform.NetTranslation().y();
    zoffset=pTransform.NetTranslation().z();
 
    G4ThreeVector pt[8];   // vertices after translation
    pt[0]=G4ThreeVector(xoffset-fDz*fTthetaCphi-fDy*fTalpha-fDx,
                        yoffset-fDz*fTthetaSphi-fDy,zoffset-fDz);
    pt[1]=G4ThreeVector(xoffset-fDz*fTthetaCphi-fDy*fTalpha+fDx,
                        yoffset-fDz*fTthetaSphi-fDy,zoffset-fDz);
    pt[2]=G4ThreeVector(xoffset-fDz*fTthetaCphi+fDy*fTalpha-fDx,
                        yoffset-fDz*fTthetaSphi+fDy,zoffset-fDz);
    pt[3]=G4ThreeVector(xoffset-fDz*fTthetaCphi+fDy*fTalpha+fDx,
                        yoffset-fDz*fTthetaSphi+fDy,zoffset-fDz);
    pt[4]=G4ThreeVector(xoffset+fDz*fTthetaCphi-fDy*fTalpha-fDx,
                        yoffset+fDz*fTthetaSphi-fDy,zoffset+fDz);
    pt[5]=G4ThreeVector(xoffset+fDz*fTthetaCphi-fDy*fTalpha+fDx,
                        yoffset+fDz*fTthetaSphi-fDy,zoffset+fDz);
    pt[6]=G4ThreeVector(xoffset+fDz*fTthetaCphi+fDy*fTalpha-fDx,
                        yoffset+fDz*fTthetaSphi+fDy,zoffset+fDz);
    pt[7]=G4ThreeVector(xoffset+fDz*fTthetaCphi+fDy*fTalpha+fDx,
                        yoffset+fDz*fTthetaSphi+fDy,zoffset+fDz);
    zMin=zoffset-fDz;
    zMax=zoffset+fDz;
    if ( pVoxelLimit.IsZLimited() )
    {
      if   ( (zMin>pVoxelLimit.GetMaxZExtent()+kCarTolerance)
          || (zMax<pVoxelLimit.GetMinZExtent()-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (zMin<pVoxelLimit.GetMinZExtent())
        {
          zMin=pVoxelLimit.GetMinZExtent();
        }
        if (zMax>pVoxelLimit.GetMaxZExtent())
        {
          zMax=pVoxelLimit.GetMaxZExtent();
        }
      }
    }

    temp[0] = pt[0].y()+(pt[4].y()-pt[0].y())
                       *(zMin-pt[0].z())/(pt[4].z()-pt[0].z()) ;
    temp[1] = pt[0].y()+(pt[4].y()-pt[0].y())
                       *(zMax-pt[0].z())/(pt[4].z()-pt[0].z()) ;
    temp[2] = pt[2].y()+(pt[6].y()-pt[2].y())
                       *(zMin-pt[2].z())/(pt[6].z()-pt[2].z()) ;
    temp[3] = pt[2].y()+(pt[6].y()-pt[2].y())
                       *(zMax-pt[2].z())/(pt[6].z()-pt[2].z()) ;        
    yMax = yoffset - std::fabs(fDz*fTthetaSphi) - fDy - fDy ;
    yMin = -yMax ;
    for(i=0;i<4;i++)
    {
      if(temp[i] > yMax) yMax = temp[i] ;
      if(temp[i] < yMin) yMin = temp[i] ;
    }
      
    if (pVoxelLimit.IsYLimited())
    {
      if ( (yMin>pVoxelLimit.GetMaxYExtent()+kCarTolerance)
        || (yMax<pVoxelLimit.GetMinYExtent()-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (yMin<pVoxelLimit.GetMinYExtent())
        {
          yMin=pVoxelLimit.GetMinYExtent();
        }
        if (yMax>pVoxelLimit.GetMaxYExtent())
        {
          yMax=pVoxelLimit.GetMaxYExtent();
        }
      }
    }

    temp[0] = pt[0].x()+(pt[4].x()-pt[0].x())
                       *(zMin-pt[0].z())/(pt[4].z()-pt[0].z()) ;
    temp[1] = pt[0].x()+(pt[4].x()-pt[0].x())
                       *(zMax-pt[0].z())/(pt[4].z()-pt[0].z()) ;
    temp[2] = pt[2].x()+(pt[6].x()-pt[2].x())
                       *(zMin-pt[2].z())/(pt[6].z()-pt[2].z()) ;
    temp[3] = pt[2].x()+(pt[6].x()-pt[2].x())
                       *(zMax-pt[2].z())/(pt[6].z()-pt[2].z()) ;
    temp[4] = pt[3].x()+(pt[7].x()-pt[3].x())
                       *(zMin-pt[3].z())/(pt[7].z()-pt[3].z()) ;
    temp[5] = pt[3].x()+(pt[7].x()-pt[3].x())
                       *(zMax-pt[3].z())/(pt[7].z()-pt[3].z()) ;
    temp[6] = pt[1].x()+(pt[5].x()-pt[1].x())
                       *(zMin-pt[1].z())/(pt[5].z()-pt[1].z()) ;
    temp[7] = pt[1].x()+(pt[5].x()-pt[1].x())
                       *(zMax-pt[1].z())/(pt[5].z()-pt[1].z()) ;

    xMax = xoffset - std::fabs(fDz*fTthetaCphi) - fDx - fDx -fDx - fDx;
    xMin = -xMax ;
    for(i=0;i<8;i++)
    {
      if(temp[i] > xMax) xMax = temp[i] ;
      if(temp[i] < xMin) xMin = temp[i] ;
    }
      // xMax/Min = f(yMax/Min) ?
    if (pVoxelLimit.IsXLimited())
    {
      if ( (xMin>pVoxelLimit.GetMaxXExtent()+kCarTolerance)
        || (xMax<pVoxelLimit.GetMinXExtent()-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (xMin<pVoxelLimit.GetMinXExtent())
        {
          xMin=pVoxelLimit.GetMinXExtent();
        }
        if (xMax>pVoxelLimit.GetMaxXExtent())
        {
          xMax=pVoxelLimit.GetMaxXExtent();
        }
      }
    }

    switch (pAxis)
    {
      case kXAxis:
        pMin=xMin;
        pMax=xMax;
        break;
      case kYAxis:
        pMin=yMin;
        pMax=yMax;
        break;
      case kZAxis:
        pMin=zMin;
        pMax=zMax;
        break;
      default:
        break;
    }

    pMin-=kCarTolerance;
    pMax+=kCarTolerance;
    flag = true;
  }
  else
  {
    // General rotated case - create and clip mesh to boundaries

    G4bool existsAfterClip=false;
    G4ThreeVectorList *vertices;

    pMin=+kInfinity;
    pMax=-kInfinity;

    // Calculate rotated vertex coordinates

    vertices=CreateRotatedVertices(pTransform);
    ClipCrossSection(vertices,0,pVoxelLimit,pAxis,pMin,pMax);
    ClipCrossSection(vertices,4,pVoxelLimit,pAxis,pMin,pMax);
    ClipBetweenSections(vertices,0,pVoxelLimit,pAxis,pMin,pMax);
      
    if (pMin!=kInfinity||pMax!=-kInfinity)
    {
      existsAfterClip=true;
        
      // Add 2*tolerance to avoid precision troubles
      //
      pMin-=kCarTolerance;
      pMax+=kCarTolerance;
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
        (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5);
        
      if (Inside(pTransform.Inverse().TransformPoint(clipCentre))!=kOutside)
      {
        existsAfterClip=true;
        pMin=pVoxelLimit.GetMinExtent(pAxis);
        pMax=pVoxelLimit.GetMaxExtent(pAxis);
      }
    }
    delete vertices ;          //  'new' in the function called
    flag = existsAfterClip ;
  }
  return flag;
}

/////////////////////////////////////////////////////////////////////////////
//
// Check in p is inside/on surface/outside solid

EInside G4Para::Inside( const G4ThreeVector& p ) const
{
  G4double xt, yt, yt1;
  EInside  in = kOutside;

  yt1 = p.y() - fTthetaSphi*p.z();
  yt  = std::fabs(yt1) ;

  // xt = std::fabs( p.x() - fTthetaCphi*p.z() - fTalpha*yt );

  xt = std::fabs( p.x() - fTthetaCphi*p.z() - fTalpha*yt1 );

  if ( std::fabs( p.z() ) <= fDz - kCarTolerance*0.5)
  {
    if (yt <= fDy - kCarTolerance*0.5)
    {
      if      ( xt <= fDx - kCarTolerance*0.5 ) in = kInside;
      else if ( xt <= fDx + kCarTolerance*0.5 ) in = kSurface;
    }
    else if ( yt <= fDy + kCarTolerance*0.5)
    {
      if ( xt <= fDx + kCarTolerance*0.5 ) in = kSurface;  
    }
  }
  else  if ( std::fabs(p.z()) <= fDz + kCarTolerance*0.5 )
  {
    if ( yt <= fDy + kCarTolerance*0.5)
    {
      if ( xt <= fDx + kCarTolerance*0.5 ) in = kSurface;  
    }
  }
  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If 2+ sides equidistant, first side's normal returned (arbitrarily)

G4ThreeVector G4Para::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4ThreeVector norm, sumnorm(0.,0.,0.);
  G4int noSurfaces = 0; 
  G4double distx,disty,distz;
  G4double newpx,newpy,xshift;
  G4double calpha,salpha;      // Sin/Cos(alpha) - needed to recalc G4Parameter
  G4double tntheta,cosntheta;  // tan and cos of normal's theta component
  G4double ycomp;
  G4double delta = 0.5*kCarTolerance;

  newpx  = p.x()-fTthetaCphi*p.z();
  newpy  = p.y()-fTthetaSphi*p.z();

  calpha = 1/std::sqrt(1+fTalpha*fTalpha);
  if (fTalpha) {salpha = -calpha/fTalpha;} // NOTE: using MINUS std::sin(alpha)
  else         {salpha = 0.;}
  
  //  xshift = newpx*calpha+newpy*salpha;
  xshift = newpx - newpy*fTalpha;

  //  distx  = std::fabs(std::fabs(xshift)-fDx*calpha);
  distx  = std::fabs(std::fabs(xshift)-fDx);
  disty  = std::fabs(std::fabs(newpy)-fDy);
  distz  = std::fabs(std::fabs(p.z())-fDz);

  tntheta   = fTthetaCphi*calpha + fTthetaSphi*salpha;
  cosntheta = 1/std::sqrt(1+tntheta*tntheta);
  ycomp     = 1/std::sqrt(1+fTthetaSphi*fTthetaSphi);

  G4ThreeVector nX  = G4ThreeVector( calpha*cosntheta,
                                     salpha*cosntheta,
                                    -tntheta*cosntheta);
  G4ThreeVector nY  = G4ThreeVector( 0, ycomp,-fTthetaSphi*ycomp);
  G4ThreeVector nZ  = G4ThreeVector( 0, 0,  1.0);

  if (distx <= delta)      
  {
    noSurfaces ++;
    if ( xshift >= 0.) {sumnorm += nX;}
    else               {sumnorm -= nX;}
  }
  if (disty <= delta)
  {
    noSurfaces ++;
    if ( newpy >= 0.)  {sumnorm += nY;}
    else               {sumnorm -= nY;}
  }
  if (distz <= delta)  
  {
    noSurfaces ++;
    if ( p.z() >= 0.)  {sumnorm += nZ;}
    else               {sumnorm -= nZ;}
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4Para::SurfaceNormal(p)", "Notification", JustWarning, 
                "Point p is not on surface !?" );
#endif 
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 ) {norm = sumnorm;}
  else                        {norm = sumnorm.unit();}

  return norm;
}


////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Para::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  ENSide  side;
  G4ThreeVector norm;
  G4double distx,disty,distz;
  G4double newpx,newpy,xshift;
  G4double calpha,salpha;  // Sin/Cos(alpha) - needed to recalc G4Parameter 
  G4double tntheta,cosntheta;  // tan and cos of normal's theta component
  G4double ycomp;

  newpx=p.x()-fTthetaCphi*p.z();
  newpy=p.y()-fTthetaSphi*p.z();

  calpha=1/std::sqrt(1+fTalpha*fTalpha);
  if (fTalpha)
  {
    salpha=-calpha/fTalpha;  // NOTE: actually use MINUS std::sin(alpha)
  }
  else
  {
    salpha=0;
  }

  xshift=newpx*calpha+newpy*salpha;

  distx=std::fabs(std::fabs(xshift)-fDx*calpha);
  disty=std::fabs(std::fabs(newpy)-fDy);
  distz=std::fabs(std::fabs(p.z())-fDz);
    
  if (distx<disty)
  {
    if (distx<distz) {side=kNX;}
    else {side=kNZ;}
  }
  else
  {
    if (disty<distz) {side=kNY;}
    else {side=kNZ;}
  }

  switch (side)
  {
    case kNX:
      tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
      if (xshift<0)
      {
        cosntheta=-1/std::sqrt(1+tntheta*tntheta);
      }
      else
      {
        cosntheta=1/std::sqrt(1+tntheta*tntheta);
      }
      norm=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
      break;
    case kNY:
      if (newpy<0)
      {
        ycomp=-1/std::sqrt(1+fTthetaSphi*fTthetaSphi);
      }
      else
      {
        ycomp=1/std::sqrt(1+fTthetaSphi*fTthetaSphi);
      }
      norm=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
      break;
    case kNZ:           // Closest to Z
      if (p.z()>=0)
      {
        norm=G4ThreeVector(0,0,1);
      }
      else
      {
        norm=G4ThreeVector(0,0,-1);
      }
      break;
  }
  return norm;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside
// - return kInfinity if no intersection
//
// ALGORITHM:
// For each component, calculate pair of minimum and maximum intersection
// values for which the particle is in the extent of the shape
// - The smallest (MAX minimum) allowed distance of the pairs is intersect
// - Z plane intersectin uses tolerance
// - XZ YZ planes use logic & *SLIGHTLY INCORRECT* tolerance
//   (this saves at least 1 sqrt, 1 multiply and 1 divide... in applicable
//    cases)
// - Note: XZ and YZ planes each divide space into four regions,
//   characterised by ss1 ss2

G4double G4Para::DistanceToIn( const G4ThreeVector& p,
                               const G4ThreeVector& v ) const
{
  G4double snxt;    // snxt = default return value
  G4double smin,smax;
  G4double tmin,tmax;
  G4double yt,vy,xt,vx;
  G4double max;
  //
  // Z Intersection range
  //
  if (v.z()>0)
  {
    max=fDz-p.z();
    if (max>kCarTolerance*0.5)
    {
      smax=max/v.z();
      smin=(-fDz-p.z())/v.z();
    }
    else
    {
      return snxt=kInfinity;
    }
  }
  else if (v.z()<0)
  {
    max=-fDz-p.z();
    if (max<-kCarTolerance*0.5)
    {
      smax=max/v.z();
      smin=(fDz-p.z())/v.z();
    }
    else
    {
      return snxt=kInfinity;
    }
  }
  else
  {
    if (std::fabs(p.z())<=fDz) // Inside
    {
      smin=0;
      smax=kInfinity;
    }
    else
    {
      return snxt=kInfinity;
    }
  }
    
  //
  // Y G4Parallel planes intersection
  //

  yt=p.y()-fTthetaSphi*p.z();
  vy=v.y()-fTthetaSphi*v.z();

  if (vy>0)
  {
    max=fDy-yt;
    if (max>kCarTolerance*0.5)
    {
      tmax=max/vy;
      tmin=(-fDy-yt)/vy;
    }
    else
    {
      return snxt=kInfinity;
    }
  }
  else if (vy<0)
  {
    max=-fDy-yt;
    if (max<-kCarTolerance*0.5)
    {
      tmax=max/vy;
      tmin=(fDy-yt)/vy;
    }
    else
    {
      return snxt=kInfinity;
    }
  }
  else
  {
    if (std::fabs(yt)<=fDy)
    {
      tmin=0;
      tmax=kInfinity;
    }
    else
    {
      return snxt=kInfinity;
    }
  }

  // Re-Calc valid intersection range
  //
  if (tmin>smin) smin=tmin;
  if (tmax<smax) smax=tmax;
  if (smax<=smin)
  {
    return snxt=kInfinity;
  }
  else
  {
    //
    // X G4Parallel planes intersection
    //
    xt=p.x()-fTthetaCphi*p.z()-fTalpha*yt;
    vx=v.x()-fTthetaCphi*v.z()-fTalpha*vy;
    if (vx>0)
    {
      max=fDx-xt;
      if (max>kCarTolerance*0.5)
      {
        tmax=max/vx;
        tmin=(-fDx-xt)/vx;
      }
      else
      {
        return snxt=kInfinity;
      }
    }
    else if (vx<0)
    {
      max=-fDx-xt;
      if (max<-kCarTolerance*0.5)
      {
        tmax=max/vx;
        tmin=(fDx-xt)/vx;
      }
      else
      {
        return snxt=kInfinity;
      }
    }
    else
    {
      if (std::fabs(xt)<=fDx)
      {
        tmin=0;
        tmax=kInfinity;
      }
      else
      {
        return snxt=kInfinity;
      }
    }
    if (tmin>smin) smin=tmin;
    if (tmax<smax) smax=tmax;
  }

  if (smax>0&&smin<smax)
  {
    if (smin>0)
    {
      snxt=smin;
    }
    else
    {
      snxt=0;
    }
  }
  else
  {
    snxt=kInfinity;
  }
  return snxt;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// - Returns 0 is point inside

G4double G4Para::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe=0.0;
  G4double distz1,distz2,disty1,disty2,distx1,distx2;
  G4double trany,cosy,tranx,cosx;

  // Z planes
  //
  distz1=p.z()-fDz;
  distz2=-fDz-p.z();
  if (distz1>distz2)
  {
    safe=distz1;
  }
  else
  {
    safe=distz2;
  }

  trany=p.y()-fTthetaSphi*p.z(); // Transformed y into `box' system

  // Transformed x into `box' system
  //
  cosy=1.0/std::sqrt(1.0+fTthetaSphi*fTthetaSphi);
  disty1=(trany-fDy)*cosy;
  disty2=(-fDy-trany)*cosy;
    
  if (disty1>safe) safe=disty1;
  if (disty2>safe) safe=disty2;

  tranx=p.x()-fTthetaCphi*p.z()-fTalpha*trany;
  cosx=1.0/std::sqrt(1.0+fTalpha*fTalpha+fTthetaCphi*fTthetaCphi);
  distx1=(tranx-fDx)*cosx;
  distx2=(-fDx-tranx)*cosx;
    
  if (distx1>safe) safe=distx1;
  if (distx2>safe) safe=distx2;
    
  if (safe<0) safe=0;
  return safe;  
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance

G4double G4Para::DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                               const G4bool calcNorm,
                               G4bool *validNorm, G4ThreeVector *n) const
{
  ESide side = kUndef;
  G4double snxt;    // snxt = return value
  G4double max,tmax;
  G4double yt,vy,xt,vx;

  G4double ycomp,calpha,salpha,tntheta,cosntheta;

  //
  // Z Intersections
  //

  if (v.z()>0)
  {
    max=fDz-p.z();
    if (max>kCarTolerance*0.5)
    {
      snxt=max/v.z();
      side=kPZ;
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(0,0,1);
      }
      return snxt=0;
    }
  }
  else if (v.z()<0)
  {
    max=-fDz-p.z();
    if (max<-kCarTolerance*0.5)
    {
      snxt=max/v.z();
      side=kMZ;
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true;
        *n=G4ThreeVector(0,0,-1);
      }
      return snxt=0;
    }
  }
  else
  {
    snxt=kInfinity;
  }
    
  //
  // Y plane intersection
  //

  yt=p.y()-fTthetaSphi*p.z();
  vy=v.y()-fTthetaSphi*v.z();

  if (vy>0)
  {
    max=fDy-yt;
    if (max>kCarTolerance*0.5)
    {
      tmax=max/vy;
      if (tmax<snxt)
      {
        snxt=tmax;
        side=kPY;
      }
    }
    else
    {
      if (calcNorm)
      {      
        *validNorm=true; // Leaving via plus Y
        ycomp=1/std::sqrt(1+fTthetaSphi*fTthetaSphi);
        *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
      }
      return snxt=0;
    }
  }
  else if (vy<0)
  {
    max=-fDy-yt;
    if (max<-kCarTolerance*0.5)
    {
      tmax=max/vy;
      if (tmax<snxt)
      {
        snxt=tmax;
        side=kMY;
      }
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true; // Leaving via minus Y
        ycomp=-1/std::sqrt(1+fTthetaSphi*fTthetaSphi);
        *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
      }
      return snxt=0;
    }
  }

  //
  // X plane intersection
  //

  xt=p.x()-fTthetaCphi*p.z()-fTalpha*yt;
  vx=v.x()-fTthetaCphi*v.z()-fTalpha*vy;
  if (vx>0)
  {
    max=fDx-xt;
    if (max>kCarTolerance*0.5)
    {
      tmax=max/vx;
      if (tmax<snxt)
      {
        snxt=tmax;
        side=kPX;
      }
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true; // Leaving via plus X
        calpha=1/std::sqrt(1+fTalpha*fTalpha);
        if (fTalpha)
        {
          salpha=-calpha/fTalpha;  // NOTE: actually use MINUS std::sin(alpha)
        }
        else
        {
          salpha=0;
        }
        tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
        cosntheta=1/std::sqrt(1+tntheta*tntheta);
        *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
      }
      return snxt=0;
    }
  }
  else if (vx<0)
  {
    max=-fDx-xt;
    if (max<-kCarTolerance*0.5)
    {
      tmax=max/vx;
      if (tmax<snxt)
      {
        snxt=tmax;
        side=kMX;
      }
    }
    else
    {
      if (calcNorm)
      {
        *validNorm=true; // Leaving via minus X
        calpha=1/std::sqrt(1+fTalpha*fTalpha);
        if (fTalpha)
        {
          salpha=-calpha/fTalpha;  // NOTE: actually use MINUS std::sin(alpha)
        }
        else
        {
          salpha=0;
        }
        tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
        cosntheta=-1/std::sqrt(1+tntheta*tntheta);
        *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
      }
      return snxt=0;
    }
  }

  if (calcNorm)
  {
    *validNorm=true;
    switch (side)
    {
      case kMZ:
        *n=G4ThreeVector(0,0,-1);
        break;
      case kPZ:
        *n=G4ThreeVector(0,0,1);
        break;
      case kMY:
        ycomp=-1/std::sqrt(1+fTthetaSphi*fTthetaSphi);
        *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
        break;        
      case kPY:
        ycomp=1/std::sqrt(1+fTthetaSphi*fTthetaSphi);
        *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
        break;        
      case kMX:
        calpha=1/std::sqrt(1+fTalpha*fTalpha);
        if (fTalpha)
        {
          salpha=-calpha/fTalpha;  // NOTE: actually use MINUS std::sin(alpha)
        }
        else
        {
          salpha=0;
        }
        tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
        cosntheta=-1/std::sqrt(1+tntheta*tntheta);
        *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
        break;
      case kPX:
        calpha=1/std::sqrt(1+fTalpha*fTalpha);
        if (fTalpha)
        {
          salpha=-calpha/fTalpha;  // NOTE: actually use MINUS std::sin(alpha)
        }
        else
        {
          salpha=0;
        }
        tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
        cosntheta=1/std::sqrt(1+tntheta*tntheta);
        *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
        break;
      default:
        DumpInfo();
        G4Exception("G4Para::DistanceToOut(p,v,..)","Notification",JustWarning,
                    "Undefined side for valid surface normal to solid.");
        break;
    }
  }
  return snxt;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is point outside

G4double G4Para::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0;
  G4double distz1,distz2,disty1,disty2,distx1,distx2;
  G4double trany,cosy,tranx,cosx;

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
     G4Exception("G4Para::DistanceToOut(p)", "Notification",
                 JustWarning, "Point p is outside !?" );
  }
#endif

  // Z planes
  //
  distz1=fDz-p.z();
  distz2=fDz+p.z();
  if (distz1<distz2)
  {
    safe=distz1;
  }
  else
  {
    safe=distz2;
  }

  trany=p.y()-fTthetaSphi*p.z(); // Transformed y into `box' system

  // Transformed x into `box' system
  //
  cosy=1.0/std::sqrt(1.0+fTthetaSphi*fTthetaSphi);
  disty1=(fDy-trany)*cosy;
  disty2=(fDy+trany)*cosy;
    
  if (disty1<safe) safe=disty1;
  if (disty2<safe) safe=disty2;

  tranx=p.x()-fTthetaCphi*p.z()-fTalpha*trany;
  cosx=1.0/std::sqrt(1.0+fTalpha*fTalpha+fTthetaCphi*fTthetaCphi);
  distx1=(fDx-tranx)*cosx;
  distx2=(fDx+tranx)*cosx;
    
  if (distx1<safe) safe=distx1;
  if (distx2<safe) safe=distx2;
    
  if (safe<0) safe=0;
  return safe;  
}

////////////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility

G4ThreeVectorList*
G4Para::CreateRotatedVertices( const G4AffineTransform& pTransform ) const
{
  G4ThreeVectorList *vertices;
  vertices=new G4ThreeVectorList();
  if (vertices)
  {
    vertices->reserve(8);
    G4ThreeVector vertex0(-fDz*fTthetaCphi-fDy*fTalpha-fDx,
                          -fDz*fTthetaSphi-fDy, -fDz);
    G4ThreeVector vertex1(-fDz*fTthetaCphi-fDy*fTalpha+fDx,
                          -fDz*fTthetaSphi-fDy, -fDz);
    G4ThreeVector vertex2(-fDz*fTthetaCphi+fDy*fTalpha-fDx,
                          -fDz*fTthetaSphi+fDy, -fDz);
    G4ThreeVector vertex3(-fDz*fTthetaCphi+fDy*fTalpha+fDx,
                          -fDz*fTthetaSphi+fDy, -fDz);
    G4ThreeVector vertex4(+fDz*fTthetaCphi-fDy*fTalpha-fDx,
                          +fDz*fTthetaSphi-fDy, +fDz);
    G4ThreeVector vertex5(+fDz*fTthetaCphi-fDy*fTalpha+fDx,
                          +fDz*fTthetaSphi-fDy, +fDz);
    G4ThreeVector vertex6(+fDz*fTthetaCphi+fDy*fTalpha-fDx,
                          +fDz*fTthetaSphi+fDy, +fDz);
    G4ThreeVector vertex7(+fDz*fTthetaCphi+fDy*fTalpha+fDx,
                          +fDz*fTthetaSphi+fDy, +fDz);

    vertices->push_back(pTransform.TransformPoint(vertex0));
    vertices->push_back(pTransform.TransformPoint(vertex1));
    vertices->push_back(pTransform.TransformPoint(vertex2));
    vertices->push_back(pTransform.TransformPoint(vertex3));
    vertices->push_back(pTransform.TransformPoint(vertex4));
    vertices->push_back(pTransform.TransformPoint(vertex5));
    vertices->push_back(pTransform.TransformPoint(vertex6));
    vertices->push_back(pTransform.TransformPoint(vertex7));
  }
  else
  {
    DumpInfo();
    G4Exception("G4Para::CreateRotatedVertices()",
                "FatalError", FatalException,
                "Error in allocation of vertices. Out of memory !");
  }
  return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Para::GetEntityType() const
{
  return G4String("G4Para");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Para::Clone() const
{
  return new G4Para(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Para::StreamInfo( std::ostream& os ) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Para\n"
     << " Parameters: \n"
     << "    half length X: " << fDx/mm << " mm \n"
     << "    half length Y: " << fDy/mm << " mm \n"
     << "    half length Z: " << fDz/mm << " mm \n"
     << "    std::tan(alpha)         : " << fTalpha/degree << " degrees \n"
     << "    std::tan(theta)*std::cos(phi): " << fTthetaCphi/degree
     << " degrees \n"
     << "    std::tan(theta)*std::sin(phi): " << fTthetaSphi/degree
     << " degrees \n"
     << "-----------------------------------------------------------\n";

  return os;
}

//////////////////////////////////////////////////////////////////////////////
//
// GetPointOnPlane
// Auxiliary method for Get Point on Surface
//

G4ThreeVector G4Para::GetPointOnPlane(G4ThreeVector p0, G4ThreeVector p1, 
                                      G4ThreeVector p2, G4ThreeVector p3, 
                                      G4double& area) const
{
  G4double lambda1, lambda2, chose, aOne, aTwo;
  G4ThreeVector t, u, v, w, Area, normal;
  
  t = p1 - p0;
  u = p2 - p1;
  v = p3 - p2;
  w = p0 - p3;

  Area = G4ThreeVector(w.y()*v.z() - w.z()*v.y(),
                       w.z()*v.x() - w.x()*v.z(),
                       w.x()*v.y() - w.y()*v.x());
  
  aOne = 0.5*Area.mag();
  
  Area = G4ThreeVector(t.y()*u.z() - t.z()*u.y(),
                       t.z()*u.x() - t.x()*u.z(),
                       t.x()*u.y() - t.y()*u.x());
  
  aTwo = 0.5*Area.mag();
  
  area = aOne + aTwo;
  
  chose = RandFlat::shoot(0.,aOne+aTwo);

  if( (chose>=0.) && (chose < aOne) )
  {
    lambda1 = RandFlat::shoot(0.,1.);
    lambda2 = RandFlat::shoot(0.,lambda1);
    return (p2+lambda1*v+lambda2*w);    
  }

  // else

  lambda1 = RandFlat::shoot(0.,1.);
  lambda2 = RandFlat::shoot(0.,lambda1);
  return (p0+lambda1*t+lambda2*u);    
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// Return a point (G4ThreeVector) randomly and uniformly
// selected on the solid surface

G4ThreeVector G4Para::GetPointOnSurface() const
{
  G4ThreeVector One, Two, Three, Four, Five, Six;
  G4ThreeVector pt[8] ;
  G4double chose, aOne, aTwo, aThree, aFour, aFive, aSix;

  pt[0] = G4ThreeVector(-fDz*fTthetaCphi-fDy*fTalpha-fDx,
                        -fDz*fTthetaSphi-fDy, -fDz);
  pt[1] = G4ThreeVector(-fDz*fTthetaCphi-fDy*fTalpha+fDx,
                        -fDz*fTthetaSphi-fDy, -fDz);
  pt[2] = G4ThreeVector(-fDz*fTthetaCphi+fDy*fTalpha-fDx,
                        -fDz*fTthetaSphi+fDy, -fDz);
  pt[3] = G4ThreeVector(-fDz*fTthetaCphi+fDy*fTalpha+fDx,
                        -fDz*fTthetaSphi+fDy, -fDz);
  pt[4] = G4ThreeVector(+fDz*fTthetaCphi-fDy*fTalpha-fDx,
                        +fDz*fTthetaSphi-fDy, +fDz);
  pt[5] = G4ThreeVector(+fDz*fTthetaCphi-fDy*fTalpha+fDx,
                        +fDz*fTthetaSphi-fDy, +fDz);
  pt[6] = G4ThreeVector(+fDz*fTthetaCphi+fDy*fTalpha-fDx,
                        +fDz*fTthetaSphi+fDy, +fDz);
  pt[7] = G4ThreeVector(+fDz*fTthetaCphi+fDy*fTalpha+fDx,
                        +fDz*fTthetaSphi+fDy, +fDz);

  // make sure we provide the points in a clockwise fashion

  One   = GetPointOnPlane(pt[0],pt[1],pt[3],pt[2], aOne);
  Two   = GetPointOnPlane(pt[4],pt[5],pt[7],pt[6], aTwo);
  Three = GetPointOnPlane(pt[6],pt[7],pt[3],pt[2], aThree);
  Four  = GetPointOnPlane(pt[4],pt[5],pt[1],pt[0], aFour); 
  Five  = GetPointOnPlane(pt[0],pt[2],pt[6],pt[4], aFive);
  Six   = GetPointOnPlane(pt[1],pt[3],pt[7],pt[5], aSix);

  chose = RandFlat::shoot(0.,aOne+aTwo+aThree+aFour+aFive+aSix);
  
  if( (chose>=0.) && (chose<aOne) )                    
    { return One; }
  else if(chose>=aOne && chose<aOne+aTwo)  
    { return Two; }
  else if(chose>=aOne+aTwo && chose<aOne+aTwo+aThree)
    { return Three; }
  else if(chose>=aOne+aTwo+aThree && chose<aOne+aTwo+aThree+aFour)
    { return Four; }
  else if(chose>=aOne+aTwo+aThree+aFour && chose<aOne+aTwo+aThree+aFour+aFive)
    { return Five; }
  return Six;
}

////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Para::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Para::CreatePolyhedron () const
{
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);
  G4double alpha = std::atan(fTalpha);
  G4double theta = std::atan(std::sqrt(fTthetaCphi*fTthetaCphi
                            +fTthetaSphi*fTthetaSphi));
    
  return new G4PolyhedronPara(fDx, fDy, fDz, alpha, theta, phi);
}

G4NURBS* G4Para::CreateNURBS () const
{
  // return new G4NURBSbox (fDx, fDy, fDz);
  return 0 ;
}
