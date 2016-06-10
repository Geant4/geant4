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
// $Id: G4Ellipsoid.cc 83572 2014-09-01 15:23:27Z gcosmo $
//
// class G4Ellipsoid
//
// Implementation for G4Ellipsoid class
//
// History:
//
// 10.11.99 G.Horton-Smith  -- first writing, based on G4Sphere class
// 25.02.05 G.Guerrieri -- Modified for future Geant4 release
//
// --------------------------------------------------------------------

#include "globals.hh"

#include "G4Ellipsoid.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"

#include "meshdefs.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

///////////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4Ellipsoid::G4Ellipsoid(const G4String& pName,
                               G4double pxSemiAxis,
                               G4double pySemiAxis,
                               G4double pzSemiAxis,
                               G4double pzBottomCut,
                               G4double pzTopCut)
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0),
    fCubicVolume(0.), fSurfaceArea(0.), zBottomCut(0.), zTopCut(0.)
{
  // note: for users that want to use the full ellipsoid it is useful
  // to include a default for the cuts 

  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();

  halfCarTolerance = kCarTolerance*0.5;
  halfRadTolerance = kRadTolerance*0.5;

  // Check Semi-Axis
  if ( (pxSemiAxis<=0.) || (pySemiAxis<=0.) || (pzSemiAxis<=0.) )
  {
     std::ostringstream message;
     message << "Invalid semi-axis - " << GetName();
     G4Exception("G4Ellipsoid::G4Ellipsoid()", "GeomSolids0002",
                 FatalErrorInArgument, message);
  }
  SetSemiAxis(pxSemiAxis, pySemiAxis, pzSemiAxis);

  if ( pzBottomCut == 0 && pzTopCut == 0 )
  {
     SetZCuts(-pzSemiAxis, pzSemiAxis);
  }
  else if ( (pzBottomCut < pzSemiAxis) && (pzTopCut > -pzSemiAxis)
         && (pzBottomCut < pzTopCut) )
  {
     SetZCuts(pzBottomCut, pzTopCut);
  }
  else
  {
     std::ostringstream message;
     message << "Invalid z-coordinate for cutting plane - " << GetName();
     G4Exception("G4Ellipsoid::G4Ellipsoid()", "GeomSolids0002",
                 FatalErrorInArgument, message);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Ellipsoid::G4Ellipsoid( __void__& a )
  : G4VSolid(a), fRebuildPolyhedron(false), fpPolyhedron(0), kRadTolerance(0.),
    halfCarTolerance(0.), halfRadTolerance(0.), fCubicVolume(0.),
    fSurfaceArea(0.), xSemiAxis(0.), ySemiAxis(0.), zSemiAxis(0.),
    semiAxisMax(0.), zBottomCut(0.), zTopCut(0.)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Ellipsoid::~G4Ellipsoid()
{
  delete fpPolyhedron; fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Ellipsoid::G4Ellipsoid(const G4Ellipsoid& rhs)
  : G4VSolid(rhs),
    fRebuildPolyhedron(false), fpPolyhedron(0),
    kRadTolerance(rhs.kRadTolerance),
    halfCarTolerance(rhs.halfCarTolerance),
    halfRadTolerance(rhs.halfRadTolerance),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    xSemiAxis(rhs.xSemiAxis), ySemiAxis(rhs.ySemiAxis),
    zSemiAxis(rhs.zSemiAxis), semiAxisMax(rhs.semiAxisMax),
    zBottomCut(rhs.zBottomCut), zTopCut(rhs.zTopCut)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Ellipsoid& G4Ellipsoid::operator = (const G4Ellipsoid& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   kRadTolerance = rhs.kRadTolerance;
   halfCarTolerance = rhs.halfCarTolerance;
   halfRadTolerance = rhs.halfRadTolerance;
   fCubicVolume = rhs.fCubicVolume; fSurfaceArea = rhs.fSurfaceArea;
   xSemiAxis = rhs.xSemiAxis; ySemiAxis = rhs.ySemiAxis;
   zSemiAxis = rhs.zSemiAxis; semiAxisMax = rhs.semiAxisMax;
   zBottomCut = rhs.zBottomCut; zTopCut = rhs.zTopCut;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = 0;

   return *this;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Ellipsoid::ComputeDimensions(G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4Ellipsoid::CalculateExtent(const EAxis pAxis,
                             const G4VoxelLimits& pVoxelLimit,
                             const G4AffineTransform& pTransform,
                                   G4double& pMin, G4double& pMax) const
{
  if (!pTransform.IsRotated())
  {
    // Special case handling for unrotated solid ellipsoid
    // Compute x/y/z mins and maxs for bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case

    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;

    G4double maxDiff,newMin,newMax;
    G4double xoff,yoff;

    xoffset=pTransform.NetTranslation().x();
    xMin=xoffset - xSemiAxis;
    xMax=xoffset + xSemiAxis;
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

    yoffset=pTransform.NetTranslation().y();
    yMin=yoffset - ySemiAxis;
    yMax=yoffset + ySemiAxis;
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

    zoffset=pTransform.NetTranslation().z();
    zMin=zoffset + (-zSemiAxis > zBottomCut ? -zSemiAxis : zBottomCut);
    zMax=zoffset + ( zSemiAxis < zTopCut ? zSemiAxis : zTopCut);
    if (pVoxelLimit.IsZLimited())
    {
      if ( (zMin>pVoxelLimit.GetMaxZExtent()+kCarTolerance)
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

    // if here, then known to cut bounding box around ellipsoid
    //
    xoff = (xoffset < xMin) ? (xMin-xoffset)
         : (xoffset > xMax) ? (xoffset-xMax) : 0.0;
    yoff = (yoffset < yMin) ? (yMin-yoffset)
         : (yoffset > yMax) ? (yoffset-yMax) : 0.0;

    // detailed calculations
    // NOTE: does not use X or Y offsets to adjust Z range,
    // and does not use Z offset to adjust X or Y range,
    // which is consistent with G4Sphere::CalculateExtent behavior
    //
    switch (pAxis)
    {
      case kXAxis:
        if (yoff==0.)
        {
          // YZ limits cross max/min x => no change
          //
          pMin=xMin;
          pMax=xMax;
        }
        else
        {
          // YZ limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs
          //
          maxDiff= 1.0-sqr(yoff/ySemiAxis);
          if (maxDiff < 0.0) { return false; }
          maxDiff= xSemiAxis * std::sqrt(maxDiff);
          newMin=xoffset-maxDiff;
          newMax=xoffset+maxDiff;
          pMin=(newMin<xMin) ? xMin : newMin;
          pMax=(newMax>xMax) ? xMax : newMax;
        }
        break;
      case kYAxis:
        if (xoff==0.)
        {
          // XZ limits cross max/min y => no change
          //
          pMin=yMin;
          pMax=yMax;
        }
        else
        {
          // XZ limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs
          //
          maxDiff= 1.0-sqr(xoff/xSemiAxis);
          if (maxDiff < 0.0) { return false; }
          maxDiff= ySemiAxis * std::sqrt(maxDiff);
          newMin=yoffset-maxDiff;
          newMax=yoffset+maxDiff;
          pMin=(newMin<yMin) ? yMin : newMin;
          pMax=(newMax>yMax) ? yMax : newMax;
        }
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
    return true;
  }
  else  // not rotated
  {
    G4int i,j,noEntries,noBetweenSections;
    G4bool existsAfterClip=false;

    // Calculate rotated vertex coordinates

    G4int noPolygonVertices=0;
    G4ThreeVectorList* vertices =
      CreateRotatedVertices(pTransform,noPolygonVertices);

    pMin=+kInfinity;
    pMax=-kInfinity;

    noEntries=vertices->size(); // noPolygonVertices*noPhiCrossSections
    noBetweenSections=noEntries-noPolygonVertices;
    
    G4ThreeVectorList ThetaPolygon;
    for (i=0;i<noEntries;i+=noPolygonVertices)
    {
      for(j=0;j<(noPolygonVertices/2)-1;j++)
      {
        ThetaPolygon.push_back((*vertices)[i+j]);  
        ThetaPolygon.push_back((*vertices)[i+j+1]);  
        ThetaPolygon.push_back((*vertices)[i+noPolygonVertices-2-j]);
        ThetaPolygon.push_back((*vertices)[i+noPolygonVertices-1-j]);
        CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
        ThetaPolygon.clear();
      }
    }
    for (i=0;i<noBetweenSections;i+=noPolygonVertices)
    {
      for(j=0;j<noPolygonVertices-1;j++)
      {
        ThetaPolygon.push_back((*vertices)[i+j]);  
        ThetaPolygon.push_back((*vertices)[i+j+1]);  
        ThetaPolygon.push_back((*vertices)[i+noPolygonVertices+j+1]);
        ThetaPolygon.push_back((*vertices)[i+noPolygonVertices+j]);
        CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
        ThetaPolygon.clear();
      }
      ThetaPolygon.push_back((*vertices)[i+noPolygonVertices-1]);
      ThetaPolygon.push_back((*vertices)[i]);
      ThetaPolygon.push_back((*vertices)[i+noPolygonVertices]);
      ThetaPolygon.push_back((*vertices)[i+2*noPolygonVertices-1]);
      CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
      ThetaPolygon.clear();
    }
    if ( (pMin!=kInfinity) || (pMax!=-kInfinity) )
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
      //
      G4ThreeVector
      clipCentre((pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
                 (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
                 (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5);

      if (Inside(pTransform.Inverse().TransformPoint(clipCentre))!=kOutside)
      {
        existsAfterClip=true;
        pMin=pVoxelLimit.GetMinExtent(pAxis);
        pMax=pVoxelLimit.GetMaxExtent(pAxis);
      }
    }
    delete vertices;
    return existsAfterClip;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
// Split into radius, phi, theta checks
// Each check modifies `in', or returns as approprate

EInside G4Ellipsoid::Inside(const G4ThreeVector& p) const
{
  G4double rad2oo,  // outside surface outer tolerance
           rad2oi;  // outside surface inner tolerance
  EInside in;

  // check this side of z cut first, because that's fast
  //
  if (p.z() < zBottomCut-halfRadTolerance) { return in=kOutside; }
  if (p.z() > zTopCut+halfRadTolerance)    { return in=kOutside; }

  rad2oo= sqr(p.x()/(xSemiAxis+halfRadTolerance))
        + sqr(p.y()/(ySemiAxis+halfRadTolerance))
        + sqr(p.z()/(zSemiAxis+halfRadTolerance));

  if (rad2oo > 1.0)  { return in=kOutside; }
    
  rad2oi= sqr(p.x()*(1.0+halfRadTolerance/xSemiAxis)/xSemiAxis)
      + sqr(p.y()*(1.0+halfRadTolerance/ySemiAxis)/ySemiAxis)
      + sqr(p.z()*(1.0+halfRadTolerance/zSemiAxis)/zSemiAxis);

  // Check radial surfaces
  //  sets `in' (already checked for rad2oo > 1.0)
  //
  if (rad2oi < 1.0)
  {
    in = ( (p.z() < zBottomCut+halfRadTolerance)
        || (p.z() > zTopCut-halfRadTolerance) ) ? kSurface : kInside;
    if ( rad2oi > 1.0-halfRadTolerance )  { in=kSurface; }
  }
  else 
  {
    in = kSurface;
  }
  return in;

}

///////////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p not protected against p=0

G4ThreeVector G4Ellipsoid::SurfaceNormal( const G4ThreeVector& p) const
{
  G4double distR, distZBottom, distZTop;

  // normal vector with special magnitude:  parallel to normal, units 1/length
  // norm*p == 1.0 if on surface, >1.0 if outside, <1.0 if inside
  //
  G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis),
                     p.y()/(ySemiAxis*ySemiAxis),
                     p.z()/(zSemiAxis*zSemiAxis));
  G4double radius = 1.0/norm.mag();

  // approximate distance to curved surface
  //
  distR = std::fabs( (p*norm - 1.0) * radius ) / 2.0;

  // Distance to z-cut plane
  //
  distZBottom = std::fabs( p.z() - zBottomCut );
  distZTop = std::fabs( p.z() - zTopCut );

  if ( (distZBottom < distR) || (distZTop < distR) )
  {
    return G4ThreeVector(0.,0.,(distZBottom < distZTop) ? -1.0 : 1.0);
  }
  return ( norm *= radius );
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//

G4double G4Ellipsoid::DistanceToIn( const G4ThreeVector& p,
                                    const G4ThreeVector& v  ) const
{
  G4double distMin = std::min(xSemiAxis,ySemiAxis);
  const G4double dRmax = 100.*std::min(distMin,zSemiAxis);
  distMin= kInfinity;

  // check to see if Z plane is relevant
  if (p.z() <= zBottomCut+halfCarTolerance)
  {
    if (v.z() <= 0.0) { return distMin; }
    G4double distZ = (zBottomCut - p.z()) / v.z();

    if ( (distZ > -halfRadTolerance) && (Inside(p+distZ*v) != kOutside) )
    {
      // early exit since can't intercept curved surface if we reach here
      if ( std::fabs(distZ) < halfRadTolerance ) { distZ=0.; }
      return distMin= distZ;
    }
  }
  if (p.z() >= zTopCut-halfCarTolerance)
  {
    if (v.z() >= 0.0) { return distMin;}
    G4double distZ = (zTopCut - p.z()) / v.z();
    if ( (distZ > -halfRadTolerance) && (Inside(p+distZ*v) != kOutside) )
    {
      // early exit since can't intercept curved surface if we reach here
      if ( std::fabs(distZ) < halfRadTolerance ) { distZ=0.; }
      return distMin= distZ;
    }
  }
  // if fZCut1 <= p.z() <= fZCut2, then must hit curved surface

  // now check curved surface intercept
  G4double A,B,C;

  A= sqr(v.x()/xSemiAxis) + sqr(v.y()/ySemiAxis) + sqr(v.z()/zSemiAxis);
  C= sqr(p.x()/xSemiAxis) + sqr(p.y()/ySemiAxis) + sqr(p.z()/zSemiAxis) - 1.0;
  B= 2.0 * ( p.x()*v.x()/(xSemiAxis*xSemiAxis)
           + p.y()*v.y()/(ySemiAxis*ySemiAxis)
           + p.z()*v.z()/(zSemiAxis*zSemiAxis) );

  C= B*B - 4.0*A*C;
  if (C > 0.0)
  {    
    G4double distR= (-B - std::sqrt(C)) / (2.0*A);
    G4double intZ = p.z()+distR*v.z();
    if ( (distR > halfRadTolerance)
      && (intZ >= zBottomCut-halfRadTolerance)
      && (intZ <= zTopCut+halfRadTolerance) )
    { 
      distMin = distR;
    }
    else if( (distR >- halfRadTolerance)
	    && (intZ >= zBottomCut-halfRadTolerance)
	    && (intZ <= zTopCut+halfRadTolerance) )
    {
      // p is on the curved surface, DistanceToIn returns 0 or kInfinity:
      // DistanceToIn returns 0, if second root is positive (means going inside)
      // If second root is negative, DistanceToIn returns kInfinity (outside)
      //
      distR = (-B + std::sqrt(C) ) / (2.0*A);
      if(distR>0.) { distMin=0.; }
    }
    else
    {
      distR= (-B + std::sqrt(C)) / (2.0*A);
      intZ = p.z()+distR*v.z();
      if ( (distR > halfRadTolerance)
        && (intZ >= zBottomCut-halfRadTolerance)
        && (intZ <= zTopCut+halfRadTolerance) )
      {
        G4ThreeVector norm=SurfaceNormal(p);
        if (norm.dot(v)<0.) { distMin = distR; }
      }
    }
    if ( (distMin!=kInfinity) && (distMin>dRmax) ) 
    {                    // Avoid rounding errors due to precision issues on
                         // 64 bits systems. Split long distances and recompute
      G4double fTerm = distMin-std::fmod(distMin,dRmax);
      distMin = fTerm + DistanceToIn(p+fTerm*v,v);
    }
  }
  
  if (std::fabs(distMin)<halfRadTolerance) { distMin=0.; }
  return distMin;
} 

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Return 0 if point inside

G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector& p) const
{
  G4double distR, distZ;

  // normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  //
  G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis),
                     p.y()/(ySemiAxis*ySemiAxis),
                     p.z()/(zSemiAxis*zSemiAxis));
  G4double radius= 1.0/norm.mag();

  // approximate distance to curved surface ( <= actual distance )
  //
  distR= (p*norm - 1.0) * radius / 2.0;

  // Distance to z-cut plane
  //
  distZ= zBottomCut - p.z();
  if (distZ < 0.0)
  {
    distZ = p.z() - zTopCut;
  }

  // Distance to closest surface from outside
  //
  if (distZ < 0.0)
  {
    return (distR < 0.0) ? 0.0 : distR;
  }
  else if (distR < 0.0)
  {
    return distZ;
  }
  else
  {
    return (distZ < distR) ? distZ : distR;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p,
                                    const G4ThreeVector& v,
                                    const G4bool calcNorm,
                                          G4bool *validNorm,
                                          G4ThreeVector *n  ) const
{
  G4double distMin;
  enum surface_e {kPlaneSurf, kCurvedSurf, kNoSurf} surface;
  
  distMin= kInfinity;
  surface= kNoSurf;

  // check to see if Z plane is relevant
  //
  if (v.z() < 0.0)
  {
    G4double distZ = (zBottomCut - p.z()) / v.z();
    if (distZ < 0.0)
    {
      distZ= 0.0;
      if (!calcNorm) {return 0.0;}
    }
    distMin= distZ;
    surface= kPlaneSurf;
  }
  if (v.z() > 0.0)
  {
    G4double distZ = (zTopCut - p.z()) / v.z();
    if (distZ < 0.0)
    {
      distZ= 0.0;
      if (!calcNorm) {return 0.0;}
    }
    distMin= distZ;
    surface= kPlaneSurf;
  }

  // normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  //
  G4ThreeVector nearnorm(p.x()/(xSemiAxis*xSemiAxis),
                         p.y()/(ySemiAxis*ySemiAxis),
                         p.z()/(zSemiAxis*zSemiAxis));
  
  // now check curved surface intercept
  //
  G4double A,B,C;
  
  A= sqr(v.x()/xSemiAxis) + sqr(v.y()/ySemiAxis) + sqr(v.z()/zSemiAxis);
  C= (p * nearnorm) - 1.0;
  B= 2.0 * (v * nearnorm);

  C= B*B - 4.0*A*C;
  if (C > 0.0)
  {
    G4double distR= (-B + std::sqrt(C) ) / (2.0*A);
    if (distR < 0.0)
    {
      distR= 0.0;
      if (!calcNorm) {return 0.0;}
    }
    if (distR < distMin)
    {
      distMin= distR;
      surface= kCurvedSurf;
    }
  }

  // set normal if requested
  //
  if (calcNorm)
  {
    if (surface == kNoSurf)
    {
      *validNorm = false;
    }
    else
    {
      *validNorm = true;
      switch (surface)
      {
        case kPlaneSurf:
          *n= G4ThreeVector(0.,0.,(v.z() > 0.0 ? 1. : -1.));
          break;
        case kCurvedSurf:
        {
          G4ThreeVector pexit= p + distMin*v;
          G4ThreeVector truenorm(pexit.x()/(xSemiAxis*xSemiAxis),
                                 pexit.y()/(ySemiAxis*ySemiAxis),
                                 pexit.z()/(zSemiAxis*zSemiAxis));
          truenorm *= 1.0/truenorm.mag();
          *n= truenorm;
        } break;
        default:           // Should never reach this case ...
          DumpInfo();
          std::ostringstream message;
          G4int oldprc = message.precision(16);
          message << "Undefined side for valid surface normal to solid."
                  << G4endl
                  << "Position:"  << G4endl
                  << "   p.x() = "   << p.x()/mm << " mm" << G4endl
                  << "   p.y() = "   << p.y()/mm << " mm" << G4endl
                  << "   p.z() = "   << p.z()/mm << " mm" << G4endl
                  << "Direction:" << G4endl << G4endl
                  << "   v.x() = "   << v.x() << G4endl
                  << "   v.y() = "   << v.y() << G4endl
                  << "   v.z() = "   << v.z() << G4endl
                  << "Proposed distance :" << G4endl
                  << "   distMin = "    << distMin/mm << " mm";
          message.precision(oldprc);
          G4Exception("G4Ellipsoid::DistanceToOut(p,v,..)",
                      "GeomSolids1002", JustWarning, message);
          break;
      }
    }
  }
   
  return distMin;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p) const
{
  G4double distR, distZ;

#ifdef G4SPECSDEBUG
  if( Inside(p) == kOutside )
  {
     DumpInfo();
     std::ostringstream message;
     G4int oldprc = message.precision(16);
     message << "Point p is outside !?" << G4endl
             << "Position:"  << G4endl
             << "   p.x() = "   << p.x()/mm << " mm" << G4endl
             << "   p.y() = "   << p.y()/mm << " mm" << G4endl
             << "   p.z() = "   << p.z()/mm << " mm";
     message.precision(oldprc) ;
     G4Exception("G4Ellipsoid::DistanceToOut(p)", "GeomSolids1002",
                 JustWarning, message);
  }
#endif

  // Normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  //
  G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis),
                     p.y()/(ySemiAxis*ySemiAxis),
                     p.z()/(zSemiAxis*zSemiAxis));

  // the following is a safe inlined "radius= min(1.0/norm.mag(),p.mag())
  //
  G4double radius= p.mag();
  G4double tmp= norm.mag();
  if ( (tmp > 0.0) && (1.0 < radius*tmp) ) {radius = 1.0/tmp;}

  // Approximate distance to curved surface ( <= actual distance )
  //
  distR = (1.0 - p*norm) * radius / 2.0;
    
  // Distance to z-cut plane
  //
  distZ = p.z() - zBottomCut;
  if (distZ < 0.0) {distZ= zTopCut - p.z();}

  // Distance to closest surface from inside
  //
  if ( (distZ < 0.0) || (distR < 0.0) )
  {
    return 0.0;
  }
  else
  {
    return (distZ < distR) ? distZ : distR;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

G4ThreeVectorList*
G4Ellipsoid::CreateRotatedVertices(const G4AffineTransform& pTransform,
                                         G4int& noPolygonVertices) const
{
  G4ThreeVectorList *vertices;
  G4ThreeVector vertex;
  G4double meshAnglePhi, meshRMaxFactor,
           crossAnglePhi, coscrossAnglePhi, sincrossAnglePhi, sAnglePhi;
  G4double meshTheta, crossTheta, startTheta;
  G4double rMaxX, rMaxY, rMaxZ, rMaxMax, rx, ry, rz;
  G4int crossSectionPhi, noPhiCrossSections, crossSectionTheta, noThetaSections;

  // Phi cross sections
  //
  noPhiCrossSections=G4int (twopi/kMeshAngleDefault)+1;  // = 9!
    
/*
  if (noPhiCrossSections<kMinMeshSections)        // <3
  {
    noPhiCrossSections=kMinMeshSections;
  }
  else if (noPhiCrossSections>kMaxMeshSections)   // >37
  {
    noPhiCrossSections=kMaxMeshSections;
  }
*/
  meshAnglePhi=twopi/(noPhiCrossSections-1);
    
  // Set start angle such that mesh will be at fRMax
  // on the x axis. Will give better extent calculations when not rotated.
    
  sAnglePhi = -meshAnglePhi*0.5;

  // Theta cross sections
    
  noThetaSections = G4int(pi/kMeshAngleDefault)+3;  //  = 7!

/*
  if (noThetaSections<kMinMeshSections)       // <3
  {
    noThetaSections=kMinMeshSections;
  }
  else if (noThetaSections>kMaxMeshSections)  // >37
  {
    noThetaSections=kMaxMeshSections;
  }
*/
  meshTheta= pi/(noThetaSections-2);
    
  // Set start angle such that mesh will be at fRMax
  // on the z axis. Will give better extent calculations when not rotated.
    
  startTheta = -meshTheta*0.5;

  meshRMaxFactor =  1.0/std::cos(0.5*
                    std::sqrt(meshAnglePhi*meshAnglePhi+meshTheta*meshTheta));
  rMaxMax= (xSemiAxis > ySemiAxis ? xSemiAxis : ySemiAxis);
  if (zSemiAxis > rMaxMax) rMaxMax= zSemiAxis;
  rMaxX= xSemiAxis + rMaxMax*(meshRMaxFactor-1.0);
  rMaxY= ySemiAxis + rMaxMax*(meshRMaxFactor-1.0);
  rMaxZ= zSemiAxis + rMaxMax*(meshRMaxFactor-1.0);
  G4double* cosCrossTheta = new G4double[noThetaSections];
  G4double* sinCrossTheta = new G4double[noThetaSections];    
  vertices=new G4ThreeVectorList(noPhiCrossSections*noThetaSections);
  if (vertices && cosCrossTheta && sinCrossTheta)
  {
    for (crossSectionTheta=0; crossSectionTheta<noThetaSections;
         crossSectionTheta++)
    {
      // Compute sine and cosine table (for historical reasons)
      //
      crossTheta=startTheta+crossSectionTheta*meshTheta;
      cosCrossTheta[crossSectionTheta]=std::cos(crossTheta);
      sinCrossTheta[crossSectionTheta]=std::sin(crossTheta);
    }
    for (crossSectionPhi=0; crossSectionPhi<noPhiCrossSections;
         crossSectionPhi++)
    {
      crossAnglePhi=sAnglePhi+crossSectionPhi*meshAnglePhi;
      coscrossAnglePhi=std::cos(crossAnglePhi);
      sincrossAnglePhi=std::sin(crossAnglePhi);
      for (crossSectionTheta=0; crossSectionTheta<noThetaSections;
           crossSectionTheta++)
      {
        // Compute coordinates of cross section at section crossSectionPhi
        //
        rx= sinCrossTheta[crossSectionTheta]*coscrossAnglePhi*rMaxX;
        ry= sinCrossTheta[crossSectionTheta]*sincrossAnglePhi*rMaxY;
        rz= cosCrossTheta[crossSectionTheta]*rMaxZ;
        if (rz < zBottomCut)
          { rz= zBottomCut; }
        if (rz > zTopCut)
          { rz= zTopCut; }
        vertex= G4ThreeVector(rx,ry,rz);
        vertices->push_back(pTransform.TransformPoint(vertex));
      }    // Theta forward     
    }    // Phi
    noPolygonVertices = noThetaSections ;
  }
  else
  {
    DumpInfo();
    G4Exception("G4Ellipsoid::CreateRotatedVertices()",
                "GeomSolids0003", FatalException,
                "Error in allocation of vertices. Out of memory !");
  }

  delete[] cosCrossTheta;
  delete[] sinCrossTheta;

  return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType

G4GeometryType G4Ellipsoid::GetEntityType() const
{
  return G4String("G4Ellipsoid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4Ellipsoid::Clone() const
{
  return new G4Ellipsoid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Ellipsoid::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Ellipsoid\n"
     << " Parameters: \n"

     << "    semi-axis x: " << xSemiAxis/mm << " mm \n"
     << "    semi-axis y: " << ySemiAxis/mm << " mm \n"
     << "    semi-axis z: " << zSemiAxis/mm << " mm \n"
     << "    max semi-axis: " << semiAxisMax/mm << " mm \n"
     << "    lower cut plane level z: " << zBottomCut/mm << " mm \n"
     << "    upper cut plane level z: " << zTopCut/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Ellipsoid::GetPointOnSurface() const
{
  G4double aTop, aBottom, aCurved, chose, xRand, yRand, zRand, phi;
  G4double cosphi, sinphi, costheta, sintheta, alpha, beta, max1, max2, max3;

  max1  = xSemiAxis > ySemiAxis ? xSemiAxis : ySemiAxis;
  max1  = max1 > zSemiAxis ? max1 : zSemiAxis;
  if (max1 == xSemiAxis)      { max2 = ySemiAxis; max3 = zSemiAxis; }
  else if (max1 == ySemiAxis) { max2 = xSemiAxis; max3 = zSemiAxis; }
  else                        { max2 = xSemiAxis; max3 = ySemiAxis; }

  phi   = RandFlat::shoot(0.,twopi);
  
  cosphi = std::cos(phi);   sinphi = std::sin(phi);
  costheta = RandFlat::shoot(zBottomCut,zTopCut)/zSemiAxis;
  sintheta = std::sqrt(1.-sqr(costheta));
  
  alpha = 1.-sqr(max2/max1); beta  = 1.-sqr(max3/max1);
  
  aTop    = pi*xSemiAxis*ySemiAxis*(1 - sqr(zTopCut/zSemiAxis));
  aBottom = pi*xSemiAxis*ySemiAxis*(1 - sqr(zBottomCut/zSemiAxis));
  
  // approximation
  // from:" http://www.citr.auckland.ac.nz/techreports/2004/CITR-TR-139.pdf"
  aCurved = 4.*pi*max1*max2*(1.-1./6.*(alpha+beta)-
                            1./120.*(3.*sqr(alpha)+2.*alpha*beta+3.*sqr(beta)));

  aCurved *= 0.5*(1.2*zTopCut/zSemiAxis - 1.2*zBottomCut/zSemiAxis);
  
  if( ( zTopCut >= zSemiAxis && zBottomCut <= -1.*zSemiAxis )
   || ( zTopCut == 0 && zBottomCut ==0 ) )
  {
    aTop = 0; aBottom = 0;
  }
  
  chose = RandFlat::shoot(0.,aTop + aBottom + aCurved); 
  
  if(chose < aCurved)
  { 
    xRand = xSemiAxis*sintheta*cosphi;
    yRand = ySemiAxis*sintheta*sinphi;
    zRand = zSemiAxis*costheta;
    return G4ThreeVector (xRand,yRand,zRand); 
  }
  else if(chose >= aCurved && chose < aCurved + aTop)
  {
    xRand = RandFlat::shoot(-1.,1.)*xSemiAxis
          * std::sqrt(1-sqr(zTopCut/zSemiAxis));
    yRand = RandFlat::shoot(-1.,1.)*ySemiAxis
          * std::sqrt(1.-sqr(zTopCut/zSemiAxis)-sqr(xRand/xSemiAxis));
    zRand = zTopCut;
    return G4ThreeVector (xRand,yRand,zRand);
  }
  else
  {
    xRand = RandFlat::shoot(-1.,1.)*xSemiAxis
          * std::sqrt(1-sqr(zBottomCut/zSemiAxis));
    yRand = RandFlat::shoot(-1.,1.)*ySemiAxis
          * std::sqrt(1.-sqr(zBottomCut/zSemiAxis)-sqr(xRand/xSemiAxis)); 
    zRand = zBottomCut;
    return G4ThreeVector (xRand,yRand,zRand);
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Ellipsoid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4VisExtent G4Ellipsoid::GetExtent() const
{
  // Define the sides of the box into which the G4Ellipsoid instance would fit.
  //
  return G4VisExtent (-semiAxisMax, semiAxisMax,
                      -semiAxisMax, semiAxisMax,
                      -semiAxisMax, semiAxisMax);
}

G4Polyhedron* G4Ellipsoid::CreatePolyhedron () const
{
  return new G4PolyhedronEllipsoid(xSemiAxis, ySemiAxis, zSemiAxis,
                                   zBottomCut, zTopCut);
}

G4Polyhedron* G4Ellipsoid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      G4AutoLock l(&polyhedronMutex);
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
      l.unlock();
    }
  return fpPolyhedron;
}
