//
// G4Ellipsoid.cc,v 1.0 2005/02/25 12:00 gguerrie Exp
//
// class G4Ellipsoid
//
// Implementation for G4Ellipsoid class
//
// History:
// 25.02.05 G.Guerrieri  -- first writing, based on G4Sphere class
//
// 09.06.05 G.Guerrieri -- Modified for Geant4 release
// 

#include <assert.h>

#include "globals.hh"

#include "G4Ellipsoid.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

// useful utility function
static inline G4double square(G4double x) { return x*x; }

// Destructor

G4Ellipsoid::~G4Ellipsoid()
{
   ;
}

// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4Ellipsoid::G4Ellipsoid(const G4String& pName,
		   G4double pxSemiAxis,
		   G4double pySemiAxis,
		   G4double pzSemiAxis,
	           G4double pzBottomCut,
		   G4double pzTopCut)  : G4CSGSolid(pName)
{

// Check Semi-Axis
    if (pxSemiAxis>0.
	&& pySemiAxis>0.
	&& pzSemiAxis>0. )
	{
	  SetSemiAxis(pxSemiAxis, pySemiAxis, pzSemiAxis);
	}
    else
	{
	    G4Exception("Error in G4Ellipsoid::G4Ellipsoid - invalid semi-axis");
	}

    if (pzBottomCut < pzSemiAxis && pzTopCut > -pzSemiAxis
	&& pzBottomCut < pzTopCut)
      {
	SetZCuts(pzBottomCut, pzTopCut);
      }
    else
      {
	G4Exception("Error in G4Ellipsoid::G4Ellipsoid - invalid z cut");
      }
}

// -------------------------------------------------------------------------------------

// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
void G4Ellipsoid::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

// --------------------------------------------------------------------------------------

// Calculate extent under transform and specified limit
G4bool G4Ellipsoid::CalculateExtent(const EAxis pAxis,
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
		if (xMin>pVoxelLimit.GetMaxXExtent()+kCarTolerance
		    ||xMax<pVoxelLimit.GetMinXExtent()-kCarTolerance)
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
		if (yMin>pVoxelLimit.GetMaxYExtent()+kCarTolerance
		    ||yMax<pVoxelLimit.GetMinYExtent()-kCarTolerance)
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
		if (zMin>pVoxelLimit.GetMaxZExtent()+kCarTolerance
		    ||zMax<pVoxelLimit.GetMinZExtent()-kCarTolerance)
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
	  xoff= (xoffset < xMin) ? (xMin-xoffset)
	    : (xoffset > xMax) ? (xoffset-xMax) : 0.0;
	  yoff= (yoffset < yMin) ? (yMin-yoffset)
	    : (yoffset > yMax) ? (yoffset-yMax) : 0.0;

// detailed calculations
// NOTE: does not use X or Y offsets to adjust Z range,
//  and does not use Z offset to adjust X or Y range,
//  which is consistent with G4Sphere::CalculateExtent behavior (in Geant4.1.0)
	  switch (pAxis)
	    {
	      case kXAxis:
		if (yoff==0.)
		  {
// YZ limits cross max/min x => no change
		      pMin=xMin;
		      pMax=xMax;
		  }
		else
		  {
// YZ limits don't cross max/min x => compute max delta x, hence new mins/maxs
		      maxDiff= 1.0-square(yoff/ySemiAxis);
		      if (maxDiff < 0.0)
			return false;
		      maxDiff= xSemiAxis * sqrt(maxDiff);
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
		      pMin=yMin;
		      pMax=yMax;
		  }
		else
		  {
// XZ limits don't cross max/min y => compute max delta y, hence new mins/maxs
		      maxDiff= 1.0-square(xoff/xSemiAxis);
		      if (maxDiff < 0.0)
			return false;
		      maxDiff= ySemiAxis * sqrt(maxDiff);
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
    else
      {
	    G4int i,j,noEntries,noBetweenSections;
	    G4bool existsAfterClip=false;

// Calculate rotated vertex coordinates
	    G4ThreeVectorList* vertices;
            G4int  noPolygonVertices ;
	    vertices=CreateRotatedVertices(pTransform,noPolygonVertices);

	    pMin=+kInfinity;
	    pMax=-kInfinity;

	    noEntries=vertices->size(); // noPolygonVertices*noPhiCrossSections
	    noBetweenSections=noEntries-noPolygonVertices;
	    
	    G4ThreeVectorList ThetaPolygon ;
	    for (i=0;i<noEntries;i+=noPolygonVertices)
		{
		   for(j=0;j<(noPolygonVertices/2)-1;j++)
		  {
		    ThetaPolygon.push_back((*vertices)[i+j]) ;		  
		    ThetaPolygon.push_back((*vertices)[i+j+1]) ;		  
		    ThetaPolygon.push_back((*vertices)[i+noPolygonVertices-2-j]);
		    ThetaPolygon.push_back((*vertices)[i+noPolygonVertices-1-j]);
	            CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
		    ThetaPolygon.clear() ;
		  }
		}
	    for (i=0;i<noBetweenSections;i+=noPolygonVertices)
		{
		   for(j=0;j<noPolygonVertices-1;j++)
		  {
		    ThetaPolygon.push_back((*vertices)[i+j]) ;		  
		    ThetaPolygon.push_back((*vertices)[i+j+1]) ;		  
		    ThetaPolygon.push_back((*vertices)[i+noPolygonVertices+j+1]);
		    ThetaPolygon.push_back((*vertices)[i+noPolygonVertices+j]);
	            CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
		    ThetaPolygon.clear() ;
		  }
		    ThetaPolygon.push_back((*vertices)[i+noPolygonVertices-1]);
		    ThetaPolygon.push_back((*vertices)[i]) ;
		    ThetaPolygon.push_back((*vertices)[i+noPolygonVertices]) ;
		    ThetaPolygon.push_back((*vertices)[i+2*noPolygonVertices-1]);
	            CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
		    ThetaPolygon.clear() ;
		}
	    
	    if (pMin!=kInfinity || pMax!=-kInfinity)
		{
		    existsAfterClip=true;
		    
// Add 2*tolerance to avoid precision troubles
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
	    delete vertices;
	    return existsAfterClip;
      }
}

// ---------------------------------------------------------------------------------------

// Return whether point inside/outside/on surface
// Split into radius, phi, theta checks
// Each check modifies `in', or returns as approprate

EInside G4Ellipsoid::Inside(const G4ThreeVector& p) const
{
  G4double rad2oo,  // outside surface outer tolerance
           rad2oi;  // outside surface inner tolerance
  EInside in;

  // check this side of z cut first, because that's fast
  if (p.z() < zBottomCut-kRadTolerance/2.0)
    return in=kOutside;
  if (p.z() > zTopCut+kRadTolerance/2.0)
    return in=kOutside;

  rad2oo= square(p.x()/(xSemiAxis+kRadTolerance/2.))
      + square(p.y()/(ySemiAxis+kRadTolerance/2.))
      + square(p.z()/(zSemiAxis+kRadTolerance/2.));

  if (rad2oo > 1.0)
    return in=kOutside;
    
  rad2oi= square(p.x()*(1.0+kRadTolerance/2./xSemiAxis)/xSemiAxis)
      + square(p.y()*(1.0+kRadTolerance/2./ySemiAxis)/ySemiAxis)
      + square(p.z()*(1.0+kRadTolerance/2./zSemiAxis)/zSemiAxis);

//
// Check radial surfaces
//  sets `in'
// (already checked for rad2oo > 1.0)
//
    if (rad2oi < 1.0)
      {
	in=  (p.z() < zBottomCut+kRadTolerance/2.0
	      || p.z() > zTopCut-kRadTolerance/2.0) ? kSurface : kInside;
      }
    else 
      {
	in=kSurface;
      }

    return in;
}

// -------------------------------------------------------------------------------------

// Return unit normal of surface closest to p
// not protected against p=0

G4ThreeVector G4Ellipsoid::SurfaceNormal( const G4ThreeVector& p) const
{
  G4double distR, distZBottom, distZTop;

//
// normal vector with special magnitude:  parallel to normal, units 1/length
// norm*p == 1.0 if on surface, >1.0 if outside, <1.0 if inside
//
    G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis), p.y()/(ySemiAxis*ySemiAxis), p.z()/(zSemiAxis*zSemiAxis));
    G4double radius= 1.0/norm.mag();

//
// approximate distance to curved surface
//
    distR= fabs( (p*norm - 1.0) * radius ) / 2.0;
    
//
// Distance to z-cut plane
//
    distZBottom = fabs( p.z() - zBottomCut );
    distZTop = fabs( p.z() - zTopCut );

    if (distZBottom < distR || distZTop < distR)
      {
	return G4ThreeVector(0.,0., (distZBottom < distZTop) ? -1.0 : 1.0);
      }
    else
      {
	return ( norm *= radius );
      }
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//

G4double G4Ellipsoid::DistanceToIn( const G4ThreeVector& p,
				    const G4ThreeVector& v  ) const
{
  G4double distMin;
  
  distMin= kInfinity;

  // check to see if Z plane is relevant
  if (p.z() < zBottomCut) {
    if (v.z() <= 0.0)
      return distMin;
    G4double distZ = (zBottomCut - p.z()) / v.z();
    if (distZ > kRadTolerance/2.0 && Inside(p+distZ*v) != kOutside )
      {
	// early exit since can't intercept curved surface if we reach here
	return distMin= distZ;
      }
  }
  if (p.z() > zTopCut) {
    if (v.z() >= 0.0)
      return distMin;
    G4double distZ = (zTopCut - p.z()) / v.z();
    if (distZ > kRadTolerance/2.0 && Inside(p+distZ*v) != kOutside )
      {
	// early exit since can't intercept curved surface if we reach here
	return distMin= distZ;
      }
  }
  // if fZCut1 <= p.z() <= fZCut2, then must hit curved surface

  // now check curved surface intercept
  G4double A,B,C;

  A= square(v.x()/xSemiAxis) + square(v.y()/ySemiAxis) + square(v.z()/zSemiAxis);
  C= square(p.x()/xSemiAxis) + square(p.y()/ySemiAxis) + square(p.z()/zSemiAxis) - 1.0;
  B= 2.0 * ( p.x()*v.x()/(xSemiAxis*xSemiAxis) + p.y()*v.y()/(ySemiAxis*ySemiAxis)
	     + p.z()*v.z()/(zSemiAxis*zSemiAxis) );

  C= B*B - 4.0*A*C;
  if (C > 0.0)
    {
      G4double distR= (-B - sqrt(C) ) / (2.0*A);
      G4double intZ= p.z()+distR*v.z();
      if (distR > kRadTolerance/2.0
	  && intZ >= zBottomCut-kRadTolerance/2.0
	  && intZ <= zTopCut+kRadTolerance/2.0)
	{
	  distMin= distR;
	}
      else
	{
	  distR= (-B + sqrt(C) ) / (2.0*A);
	  intZ= p.z()+distR*v.z();
	  if (distR > kRadTolerance/2.0
	      && intZ >= zBottomCut-kRadTolerance/2.0
	      && intZ <= zTopCut+kRadTolerance/2.0)
	    {
	      distMin= distR;
	    }
	}
    }

  return distMin;
} 

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Return 0 if point inside

G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector& p) const
{
  G4double distR, distZ;

//
// normal vector:  parallel to normal, magnitude 1/(characteristic radius)
//
    G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis), p.y()/(ySemiAxis*ySemiAxis), p.z()/(zSemiAxis*zSemiAxis));
    G4double radius= 1.0/norm.mag();

//
// approximate distance to curved surface ( <= actual distance )
//
    distR= (p*norm - 1.0) * radius / 2.0;
    
//
// Distance to z-cut plane
//
    distZ= zBottomCut - p.z();
    if (distZ < 0.0)
      distZ= p.z() - zTopCut;

//
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

// -------------------------------------------------------------------------------------

// Calculate distance to surface of shape from `inside', allowing for tolerance

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p,
				 const G4ThreeVector& v,
			         const G4bool calcNorm,
			         G4bool *validNorm,
				 G4ThreeVector *n       ) const
{
  G4double distMin;
  enum surface_e {kPlaneSurf, kCurvedSurf, kNoSurf} surface;
  
  distMin= kInfinity;
  surface= kNoSurf;

  // check to see if Z plane is relevant
  if (v.z() < 0.0) {
    G4double distZ = (zBottomCut - p.z()) / v.z();
    if (distZ < 0.0)
      {
	distZ= 0.0;
	if (!calcNorm)
	  return 0.0;
      }
    distMin= distZ;
    surface= kPlaneSurf;
  }
  if (v.z() > 0.0) {
    G4double distZ = (zTopCut - p.z()) / v.z();
    if (distZ < 0.0)
      {
	distZ= 0.0;
	if (!calcNorm)
	  return 0.0;
      }
    distMin= distZ;
    surface= kPlaneSurf;
  }

  // normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  G4ThreeVector nearnorm(p.x()/(xSemiAxis*xSemiAxis), p.y()/(ySemiAxis*ySemiAxis), p.z()/(zSemiAxis*zSemiAxis));
  
  // now check curved surface intercept
  G4double A,B,C;
  
  A= square(v.x()/xSemiAxis) + square(v.y()/ySemiAxis) + square(v.z()/zSemiAxis);
  C= (p * nearnorm) - 1.0;
  B= 2.0 * (v * nearnorm);

  C= B*B - 4.0*A*C;
  if (C > 0.0)
    {
      G4double distR= (-B + sqrt(C) ) / (2.0*A);
      if (distR < 0.0)
	{
	  distR= 0.0;
	  if (!calcNorm)
	    return 0.0;
	}
      if (distR < distMin)
	{
	  distMin= distR;
	  surface= kCurvedSurf;
	}
    }

  // set normal if requested
  if (calcNorm)
    {
      if (surface == kNoSurf)
	{
	  *validNorm= FALSE;
	}
      else
	{
	  *validNorm= TRUE;
	  switch (surface)
	    {
	    case kPlaneSurf:
	      *n= G4ThreeVector(0.,0.,(v.z() > 1.0 ? 1. : -1.));
	      break;
	    case kCurvedSurf:
	      {
		G4ThreeVector pexit= p + distMin*v;
		G4ThreeVector truenorm(pexit.x()/(xSemiAxis*xSemiAxis),
				       pexit.y()/(ySemiAxis*ySemiAxis),
				       pexit.z()/(zSemiAxis*zSemiAxis));
		truenorm*= 1.0/truenorm.mag();
		*n= truenorm;
	      }
	      break;
	    default:
	      G4Exception("Logic error in G4Ellipsoid::DistanceToOut!");
	      break;
	    }
	}
    }
  
  return distMin;
}

// ------------------------------------------------------------------------------------

// Calcluate distance (<=actual) to closest surface of shape from inside

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p) const
{
  G4double distR, distZ;

//
// normal vector:  parallel to normal, magnitude 1/(characteristic radius)
//
    G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis), p.y()/(ySemiAxis*ySemiAxis), p.z()/(zSemiAxis*zSemiAxis));
    // the following is a safe inlined "radius= min(1.0/norm.mag(),p.mag())
    G4double radius= p.mag();
    {
      G4double tmp= norm.mag();
      if (tmp > 0.0 && 1.0 < radius*tmp) radius= 1.0/tmp;
    }

//
// approximate distance to curved surface ( <= actual distance )
//
    distR= (1.0 - p*norm) * radius / 2.0;
    
//
// Distance to z-cut plane
//
    distZ= p.z() - zBottomCut;
    if (distZ < 0.0)
      distZ= zTopCut - p.z();

//
// Distance to closest surface from inside
//
    if (distZ < 0.0 || distR < 0.0)
      {
	return 0.0;
      }
    else
      {
	return (distZ < distR) ? distZ : distR;
      }
}

// --------------------------------------------------------------------------------------

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
  G4double meshAnglePhi,meshRMaxFactor,
    crossAnglePhi,coscrossAnglePhi,sincrossAnglePhi,sAnglePhi;
  G4double meshTheta,crossTheta,startTheta;
  G4double rMaxX,rMaxY,rMaxZ,rMaxMax, rx,ry,rz;
  G4int crossSectionPhi,noPhiCrossSections,crossSectionTheta,noThetaSections;

  // Phi cross sections
    
  // noPhiCrossSections=G4int (M_PI/kMeshAngleDefault)+1;
  noPhiCrossSections=G4int (2*M_PI/kMeshAngleDefault)+1;
    
  if (noPhiCrossSections<kMinMeshSections)
    {
      noPhiCrossSections=kMinMeshSections;
    }
  else if (noPhiCrossSections>kMaxMeshSections)
    {
      noPhiCrossSections=kMaxMeshSections;
    }
  // meshAnglePhi=M_PI/(noPhiCrossSections-1);
  meshAnglePhi=2.0*M_PI/(noPhiCrossSections-1);
    
// Set start angle such that mesh will be at fRMax
// on the x axis. Will give better extent calculations when not rotated.
    
  sAnglePhi = -meshAnglePhi*0.5;

  // Theta cross sections
    
  noThetaSections = G4int(M_PI/kMeshAngleDefault)+3;
    
  if (noThetaSections<kMinMeshSections)
    {
      noThetaSections=kMinMeshSections;
    }
  else if (noThetaSections>kMaxMeshSections)
    {
      noThetaSections=kMaxMeshSections;
    }
  meshTheta= M_PI/(noThetaSections-2);
    
// Set start angle such that mesh will be at fRMax
// on the z axis. Will give better extent calculations when not rotated.
    
  startTheta = -meshTheta*0.5;

  meshRMaxFactor =  1.0/cos(0.5*hypot(meshAnglePhi,meshTheta));
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
      for (crossSectionTheta=0;crossSectionTheta<noThetaSections;crossSectionTheta++)
	{
	  // Compute sine and cosine table (for historical reasons)
	  crossTheta=startTheta+crossSectionTheta*meshTheta;
	  cosCrossTheta[crossSectionTheta]=cos(crossTheta);
	  sinCrossTheta[crossSectionTheta]=sin(crossTheta);
	}
      for (crossSectionPhi=0;crossSectionPhi<noPhiCrossSections;crossSectionPhi++)
	{
	  crossAnglePhi=sAnglePhi+crossSectionPhi*meshAnglePhi;
	  coscrossAnglePhi=cos(crossAnglePhi);
	  sincrossAnglePhi=sin(crossAnglePhi);
	  for (crossSectionTheta=0;crossSectionTheta<noThetaSections;crossSectionTheta++)
	    {
	      // Compute coordinates of cross section at section crossSectionPhi
	      rx= sinCrossTheta[crossSectionTheta]*coscrossAnglePhi*rMaxX;
	      ry= sinCrossTheta[crossSectionTheta]*sincrossAnglePhi*rMaxY;
	      rz= cosCrossTheta[crossSectionTheta]*rMaxZ;
	      if (rz < zBottomCut)
		rz= zBottomCut;
	      if (rz > zTopCut)
		rz= zTopCut;
	      vertex= G4ThreeVector(rx,ry,rz);
	      vertices->push_back(pTransform.TransformPoint(vertex));
	    }    // Theta forward     
	}       // Phi
      noPolygonVertices = noThetaSections ;
    }
  else
    {
      G4Exception("G4Ellipsoid::CreateRotatedVertices Out of memory - Cannot alloc vertices");
    }

  delete[] cosCrossTheta;
  delete[] sinCrossTheta;

  return vertices;
}

// ---------------------------------------------------------------------------------------

void G4Ellipsoid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
 // scene.AddThis (*this);
 scene.AddSolid (*this);
}

G4VisExtent G4Ellipsoid::GetExtent() const
{
  // Define the sides of the box into which the G4Ellipsoid instance would fit.
  return G4VisExtent (-semiAxisMax, semiAxisMax, -semiAxisMax, semiAxisMax, -semiAxisMax,semiAxisMax);
}

G4NURBS* G4Ellipsoid::CreateNURBS () const
{
    return new G4NURBSbox (semiAxisMax, semiAxisMax, semiAxisMax);       // Box for now!!!
}

G4Polyhedron* G4Ellipsoid::CreatePolyhedron () const
{
    return new G4PolyhedronEllipsoid (xSemiAxis, ySemiAxis, zSemiAxis, zBottomCut, zTopCut);
}

// ******************************  End of G4Ellipsoid.cc  ***************************** //
