// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Para.cc,v 1.1 1999-01-07 16:07:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Para
//
// Implementation for G4Para class
//
// History:
// 21.03.95 P.Kent Modified for `tolerant' geom
// 31.10.96 V.Grichine Modifications according G4Box/Tubs before to commit

#include <math.h>
#include "G4Para.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

// Private enum: Not for external use   	
enum ESide {kPX,kMX,kPY,kMY,kPZ,kMZ};

// used internally for normal routine
enum ENSide {kNZ,kNX,kNY};

// Constructor - check and set half-widths

void G4Para::SetAllParameters(G4double pDx, G4double pDy, G4double pDz, 
            G4double pAlpha, G4double pTheta, G4double pPhi)
{
    if (pDx>0&&pDy>0&&pDz>0)
	{
	   fDx=pDx;
	   fDy=pDy;
	   fDz=pDz;
	   fTalpha=tan(pAlpha);
	   fTthetaCphi=tan(pTheta)*cos(pPhi);
	   fTthetaSphi=tan(pTheta)*sin(pPhi);
	}
    else
	{
	   G4Exception("Error in G4Para::SetAllParameters - Invalid Length Parameters");
	}
}

G4Para::G4Para(const G4String& pName,G4double pDx, G4double pDy, G4double pDz,
           G4double pAlpha, G4double pTheta,G4double pPhi) : G4CSGSolid(pName)
{
    if (pDx>0&&pDy>0&&pDz>0)
        {
           SetAllParameters( pDx, pDy, pDz, pAlpha, pTheta, pPhi);
	}
    else
	{
	   G4Exception("Error in G4Para::G4Para - Invalid Length Parameters");
	}
}
// ----------------------------------------------------------------------------

// Constructor - Design of trapezoid based on 8 G4ThreeVector parameters, which are its vertices.
// Checking of planarity with preparation of fPlanes[] and than calculation of other members


G4Para::G4Para( const G4String& pName,
	        const G4ThreeVector pt[8]) : G4CSGSolid(pName)
{
    if (   pt[0].z()<0 && pt[0].z()==pt[1].z() && pt[0].z()==pt[2].z() && pt[0].z()==pt[3].z()
	&& pt[4].z()>0 && pt[4].z()==pt[5].z() && pt[4].z()==pt[6].z() && pt[4].z()==pt[7].z()
	&& (pt[0].z()+pt[4].z())== 0
	&& pt[0].y()==pt[1].y() && pt[2].y()==pt[3].y()
	&& pt[4].y()==pt[5].y() && pt[6].y()==pt[7].y()
	&& (pt[0].y()+pt[2].y()+pt[4].y()+pt[6].y())==0      )
     {
	   fDz = (pt[7]).z() ;
	    
	   fDy = ((pt[2]).y()-(pt[1]).y())*0.5 ;
	   fDx = ((pt[1]).x()-(pt[0]).x())*0.5 ;
	   fDx = ((pt[3]).x()-(pt[2]).x())*0.5 ;
	   fTalpha = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/fDy ;

	   //	   fDy = ((pt[6]).y()-(pt[5]).y())*0.5 ;
	   //	   fDx = ((pt[5]).x()-(pt[4]).x())*0.5 ;
	   //	   fDx = ((pt[7]).x()-(pt[6]).x())*0.5 ;
	   //	   fTalpha = ((pt[6]).x()+(pt[7]).x()-(pt[5]).x()-(pt[4]).x())*0.25/fDy ;

	   fTthetaCphi = ((pt[4]).x()+fDy*fTalpha+fDx)/fDz ;
	   fTthetaSphi = ((pt[4]).y()+fDy)/fDz ;

	}
    else
	{
	    G4Exception("Error in G4Para::G4Para - Invalid vertice coordinates");
	}

    
}

// --------------------------------------------------------------------------------------

G4Para::~G4Para()
{
   ;
}

// -----------------------------------------------------------------------------


// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

 void G4Para::ComputeDimensions(G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}


// ----------------------------------------------------------------------

// Calculate extent under transform and specified limit

G4bool G4Para::CalculateExtent(const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                               G4double& pMin, G4double& pMax) const
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
	  G4double temp[8] ;          // some points for intersection with zMin/zMax
	  
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
	    if (pVoxelLimit.IsZLimited())
		{
		    if (zMin>pVoxelLimit.GetMaxZExtent()
			||zMax<pVoxelLimit.GetMinZExtent())
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

            temp[0] = pt[0].y()+(pt[4].y()-pt[0].y())*(zMin-pt[0].z())/(pt[4].z()-pt[0].z()) ;
       	    temp[1] = pt[0].y()+(pt[4].y()-pt[0].y())*(zMax-pt[0].z())/(pt[4].z()-pt[0].z()) ;
	    temp[2] = pt[2].y()+(pt[6].y()-pt[2].y())*(zMin-pt[2].z())/(pt[6].z()-pt[2].z()) ;
	    temp[3] = pt[2].y()+(pt[6].y()-pt[2].y())*(zMax-pt[2].z())/(pt[6].z()-pt[2].z()) ;	      
	    yMax = yoffset - fabs(fDz*fTthetaSphi) - fDy - fDy ;
	    yMin = -yMax ;
	    for(i=0;i<4;i++)
	    {
	       if(temp[i] > yMax) yMax = temp[i] ;
	       if(temp[i] < yMin) yMin = temp[i] ;
	    }
	    
	    if (pVoxelLimit.IsYLimited())
		{
		    if (yMin>pVoxelLimit.GetMaxYExtent()
			||yMax<pVoxelLimit.GetMinYExtent())
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

            temp[0] = pt[0].x()+(pt[4].x()-pt[0].x())*(zMin-pt[0].z())/(pt[4].z()-pt[0].z()) ;
       	    temp[1] = pt[0].x()+(pt[4].x()-pt[0].x())*(zMax-pt[0].z())/(pt[4].z()-pt[0].z()) ;
	    temp[2] = pt[2].x()+(pt[6].x()-pt[2].x())*(zMin-pt[2].z())/(pt[6].z()-pt[2].z()) ;
	    temp[3] = pt[2].x()+(pt[6].x()-pt[2].x())*(zMax-pt[2].z())/(pt[6].z()-pt[2].z()) ;
            temp[4] = pt[3].x()+(pt[7].x()-pt[3].x())*(zMin-pt[3].z())/(pt[7].z()-pt[3].z()) ;
       	    temp[5] = pt[3].x()+(pt[7].x()-pt[3].x())*(zMax-pt[3].z())/(pt[7].z()-pt[3].z()) ;
	    temp[6] = pt[1].x()+(pt[5].x()-pt[1].x())*(zMin-pt[1].z())/(pt[5].z()-pt[1].z()) ;
	    temp[7] = pt[1].x()+(pt[5].x()-pt[1].x())*(zMax-pt[1].z())/(pt[5].z()-pt[1].z()) ;
	    
	    xMax = xoffset - fabs(fDz*fTthetaCphi) - fDx - fDx -fDx - fDx;
	    xMin = -xMax ;
	    for(i=0;i<8;i++)
	    {
	       if(temp[i] > xMax) xMax = temp[i] ;
	       if(temp[i] < xMin) xMin = temp[i] ;
	    }
	                                          // xMax/Min = f(yMax/Min) ?
	    if (pVoxelLimit.IsXLimited())
		{
		    if (xMin>pVoxelLimit.GetMaxXExtent()
			||xMax<pVoxelLimit.GetMinXExtent())
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

// --------------------------------------------------------------------------------

// Check in p is inside/on surface/outside solid
EInside G4Para::Inside(const G4ThreeVector& p) const
{
    EInside in=kOutside;
    G4double xt,yt;
// Check
    if (fabs(p.z())<=fDz-kCarTolerance/2)
	{
	    yt=fabs(p.y()-fTthetaSphi*p.z());
	    if (yt<=fDy-kCarTolerance/2)
		{
		    xt=fabs(p.x()-fTthetaCphi*p.z()-fTalpha*yt);
		    if (xt<=fDx-kCarTolerance/2)
			{
			    in=kInside;
			}
		    else if (xt<=fDx+kCarTolerance/2)
			{
			    in=kSurface;
			}
		}
	    else if (yt<=fDy+kCarTolerance/2)
		{
		    xt=fabs(p.x()-fTthetaCphi*p.z()-fTalpha*yt);
		    if (xt<=fDx+kCarTolerance/2)
			{
			    in=kSurface;
			}
		}
	}
    else  if (fabs(p.z())<=fDz+kCarTolerance/2)
	{
	    yt=fabs(p.y()-fTthetaSphi*p.z());
	    if (yt<=fDy+kCarTolerance/2)
		{
		    xt=fabs(p.x()-fTthetaCphi*p.z()-fTalpha*yt);
		    if (xt<=fDx+kCarTolerance/2)
			{
			    in=kSurface;
			}
		}
	}
    return in;
}

// -----------------------------------------------------------------------

// Calculate side nearest to p, and return normal
// If 2+ sides equidistant, first side's normal returned (arbitrarily)

G4ThreeVector G4Para::SurfaceNormal( const G4ThreeVector& p) const
{
    ENSide  side;
    G4ThreeVector norm;
    G4double distx,disty,distz;
    G4double newpx,newpy,xshift;
    G4double calpha,salpha;	// Sin/Cos(alpha) - needed to recalc G4Parameter 
    G4double tntheta,cosntheta;	// tan and cos of normal's theta component
    G4double ycomp;

    newpx=p.x()-fTthetaCphi*p.z();
    newpy=p.y()-fTthetaSphi*p.z();

    calpha=1/sqrt(1+fTalpha*fTalpha);
    if (fTalpha)
	{
	    salpha=-calpha/fTalpha;	// NOTE: actually use MINUS sin(alpha)
	}
    else
	{
	    salpha=0;
	}

    xshift=newpx*calpha+newpy*salpha;

    distx=fabs(fabs(xshift)-fDx*calpha);
    disty=fabs(fabs(newpy)-fDy);
    distz=fabs(fabs(p.z())-fDz);
    
    if (distx<disty)
	{
	    if (distx<distz) side=kNX;
	    else side=kNZ;
	}
    else
	{
	    if (disty<distz) side=kNY;
	    else side=kNZ;
	}

    switch (side)
	{
	case kNX:
	    tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
	    if (xshift<0)
		{
		    cosntheta=-1/sqrt(1+tntheta*tntheta);
		}
	    else
		{
		    cosntheta=1/sqrt(1+tntheta*tntheta);
		}
	    norm=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
	    break;
	case kNY:
	    if (newpy<0)
		{
		    ycomp=-1/sqrt(1+fTthetaSphi*fTthetaSphi);
		}
	    else
		{
		    ycomp=1/sqrt(1+fTthetaSphi*fTthetaSphi);
		}
	    norm=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
	    break;
	case kNZ:
// Closest to Z
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

// ---------------------------------------------------------------------------

// Calculate distance to shape from outside - return kInfinity if no intersection
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

G4double G4Para::DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const
{
    G4double snxt;		// snxt = default return value
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
	    if (max>kCarTolerance/2)
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
	    if (max<-kCarTolerance/2)
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
	    if (fabs(p.z())<=fDz) // Inside
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
	    if (max>kCarTolerance/2)
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
	    if (max<-kCarTolerance/2)
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
	    if (fabs(yt)<=fDy)
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
		    if (max>kCarTolerance/2)
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
		    if (max<-kCarTolerance/2)
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
		    if (fabs(xt)<=fDx)
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

// ----------------------------------------------------------------------

// Calculate exact shortest distance to any boundary from outside
// - Returns 0 is point inside

G4double G4Para::DistanceToIn(const G4ThreeVector& p) const
{
    G4double safe;
    G4double distz1,distz2,disty1,disty2,distx1,distx2;
    G4double trany,cosy,tranx,cosx;

// Z planes
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
    cosy=1.0/sqrt(1.0+fTthetaSphi*fTthetaSphi);
    disty1=(trany-fDy)*cosy;
    disty2=(-fDy-trany)*cosy;
    
    if (disty1>safe) safe=disty1;
    if (disty2>safe) safe=disty2;

    tranx=p.x()-fTthetaCphi*p.z()-fTalpha*trany;
    cosx=1.0/sqrt(1.0+fTalpha*fTalpha+fTthetaCphi*fTthetaCphi);
    distx1=(tranx-fDx)*cosx;
    distx2=(-fDx-tranx)*cosx;
    
    if (distx1>safe) safe=distx1;
    if (distx2>safe) safe=distx2;
    
    if (safe<0) safe=0;
    return safe;	
}

// ---------------------------------------------------------------------------

// Calcluate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance

G4double G4Para::DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			     const G4bool calcNorm,
			     G4bool *validNorm,G4ThreeVector *n) const
{
    ESide side;
    G4double snxt;		// snxt = return value
    G4double max,tmax;
    G4double yt,vy,xt,vx;

    G4double ycomp,calpha,salpha,tntheta,cosntheta;
//
// Z Intersections
//
    if (v.z()>0)
	{
	    max=fDz-p.z();
	    if (max>kCarTolerance/2)
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
	    if (max<-kCarTolerance/2)
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
	    if (max>kCarTolerance/2)
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
			    ycomp=1/sqrt(1+fTthetaSphi*fTthetaSphi);
			    *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
			}
		    return snxt=0;
		}
	}
    else if (vy<0)
	{
	    max=-fDy-yt;
	    if (max<-kCarTolerance/2)
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
			    ycomp=-1/sqrt(1+fTthetaSphi*fTthetaSphi);
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
	    if (max>kCarTolerance/2)
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
			    calpha=1/sqrt(1+fTalpha*fTalpha);
			    if (fTalpha)
				{
				    salpha=-calpha/fTalpha;	// NOTE: actually use MINUS sin(alpha)
				}
			    else
				{
				    salpha=0;
				}
			    tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
			    cosntheta=1/sqrt(1+tntheta*tntheta);
			    *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
			}
		    return snxt=0;
		}
	}
    else if (vx<0)
	{
	    max=-fDx-xt;
	    if (max<-kCarTolerance/2)
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
			    calpha=1/sqrt(1+fTalpha*fTalpha);
			    if (fTalpha)
				{
				    salpha=-calpha/fTalpha;	// NOTE: actually use MINUS sin(alpha)
				}
			    else
				{
				    salpha=0;
				}
			    tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
			    cosntheta=-1/sqrt(1+tntheta*tntheta);
			    *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);		       
			    return snxt=0;
			}
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
		    ycomp=-1/sqrt(1+fTthetaSphi*fTthetaSphi);
		    *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
		    break;		    
		case kPY:
		    ycomp=1/sqrt(1+fTthetaSphi*fTthetaSphi);
		    *n=G4ThreeVector(0,ycomp,-fTthetaSphi*ycomp);
		    break;		    
		case kMX:
		    calpha=1/sqrt(1+fTalpha*fTalpha);
		    if (fTalpha)
			{
			    salpha=-calpha/fTalpha;	// NOTE: actually use MINUS sin(alpha)
			}
		    else
			{
			    salpha=0;
			}
		    tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
		    cosntheta=-1/sqrt(1+tntheta*tntheta);
		    *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
		    break;
		case kPX:
		    calpha=1/sqrt(1+fTalpha*fTalpha);
		    if (fTalpha)
			{
			    salpha=-calpha/fTalpha;	// NOTE: actually use MINUS sin(alpha)
			}
		    else
			{
			    salpha=0;
			}
		    tntheta=fTthetaCphi*calpha+fTthetaSphi*salpha;
		    cosntheta=1/sqrt(1+tntheta*tntheta);
		    *n=G4ThreeVector(calpha*cosntheta,salpha*cosntheta,-tntheta*cosntheta);
		    break;
		default:
		    G4Exception("Invalid enum in G4Para::DistanceToOut");
		    break;
		    
		}
	}
    return snxt;
}

// -------------------------------------------------------------------------

// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is point outside

G4double G4Para::DistanceToOut(const G4ThreeVector& p) const
{
    G4double safe;
    G4double distz1,distz2,disty1,disty2,distx1,distx2;
    G4double trany,cosy,tranx,cosx;

// Z planes
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
    cosy=1.0/sqrt(1.0+fTthetaSphi*fTthetaSphi);
    disty1=(fDy-trany)*cosy;
    disty2=(fDy+trany)*cosy;
    
    if (disty1<safe) safe=disty1;
    if (disty2<safe) safe=disty2;

    tranx=p.x()-fTthetaCphi*p.z()-fTalpha*trany;
    cosx=1.0/sqrt(1.0+fTalpha*fTalpha+fTthetaCphi*fTthetaCphi);
    distx1=(fDx-tranx)*cosx;
    distx2=(fDx+tranx)*cosx;
    
    if (distx1<safe) safe=distx1;
    if (distx2<safe) safe=distx2;
    
    if (safe<0) safe=0;
    return safe;	
}

// ----------------------------------------------------------------------------

// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility

G4ThreeVectorList*
   G4Para::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
    G4ThreeVectorList *vertices;
    vertices=new G4ThreeVectorList(8);
    if (vertices)
	{
	    G4ThreeVector vertex0(-fDz*fTthetaCphi-fDy*fTalpha-fDx,-fDz*fTthetaSphi-fDy,-fDz);
	    G4ThreeVector vertex1(-fDz*fTthetaCphi-fDy*fTalpha+fDx,-fDz*fTthetaSphi-fDy,-fDz);
	    G4ThreeVector vertex2(-fDz*fTthetaCphi+fDy*fTalpha-fDx,-fDz*fTthetaSphi+fDy,-fDz);
	    G4ThreeVector vertex3(-fDz*fTthetaCphi+fDy*fTalpha+fDx,-fDz*fTthetaSphi+fDy,-fDz);
	    G4ThreeVector vertex4(+fDz*fTthetaCphi-fDy*fTalpha-fDx,+fDz*fTthetaSphi-fDy,+fDz);
	    G4ThreeVector vertex5(+fDz*fTthetaCphi-fDy*fTalpha+fDx,+fDz*fTthetaSphi-fDy,+fDz);
	    G4ThreeVector vertex6(+fDz*fTthetaCphi+fDy*fTalpha-fDx,+fDz*fTthetaSphi+fDy,+fDz);
	    G4ThreeVector vertex7(+fDz*fTthetaCphi+fDy*fTalpha+fDx,+fDz*fTthetaSphi+fDy,+fDz);

	    vertices->insert(pTransform.TransformPoint(vertex0));
	    vertices->insert(pTransform.TransformPoint(vertex1));
	    vertices->insert(pTransform.TransformPoint(vertex2));
	    vertices->insert(pTransform.TransformPoint(vertex3));
	    vertices->insert(pTransform.TransformPoint(vertex4));
	    vertices->insert(pTransform.TransformPoint(vertex5));
	    vertices->insert(pTransform.TransformPoint(vertex6));
	    vertices->insert(pTransform.TransformPoint(vertex7));
	}
    else
	{
	    G4Exception("G4Para::CreateRotatedVertices Out of memory - Cannot alloc vertices");
	}
    return vertices;
}


void G4Para::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
   scene.AddThis (*this);
   //  scene.AddThis (*this);
}

G4VisExtent G4Para::GetExtent() const
{
  G4int i ;
  G4double xMin,xMax,yMin,yMax ; 
  G4ThreeVector pt[8] ;
	   pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy*fTalpha-fDx,-fDz*fTthetaSphi-fDy,-fDz);
           pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy*fTalpha+fDx,-fDz*fTthetaSphi-fDy,-fDz);
	   pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy*fTalpha-fDx,-fDz*fTthetaSphi+fDy,-fDz);
	   pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy*fTalpha+fDx,-fDz*fTthetaSphi+fDy,-fDz);
	   pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy*fTalpha-fDx,+fDz*fTthetaSphi-fDy,+fDz);
	   pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy*fTalpha+fDx,+fDz*fTthetaSphi-fDy,+fDz);
	   pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy*fTalpha-fDx,+fDz*fTthetaSphi+fDy,+fDz);
	   pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy*fTalpha+fDx,+fDz*fTthetaSphi+fDy,+fDz);
	   xMax = -fabs(fDz*fTthetaCphi)-fDx-fDx-fDx-fDx ;
	   xMin = -xMax ;
	   yMax = -fabs(fDz*fTthetaSphi)-fDy-fDy ;
           yMin = -yMax ;
	   for(i=0;i<8;i++)
	   {
	      if(pt[i].x() > xMax) xMax = pt[i].x() ;
	      if(pt[i].x() < xMin) xMin = pt[i].x() ;
	      if(pt[i].y() > yMax) yMax = pt[i].y() ;
	      if(pt[i].y() < yMin) yMin = pt[i].y() ;
	   }	   
  return G4VisExtent (xMin, xMax, yMin, yMax, -fDz, fDz);
}

G4Polyhedron* G4Para::CreatePolyhedron () const
{
    G4double phi = atan2(fTthetaSphi, fTthetaCphi);
    G4double alpha = atan(fTalpha);
    G4double theta = atan(sqrt(fTthetaCphi*fTthetaCphi
                               +fTthetaSphi*fTthetaSphi));
    
    return new G4PolyhedronPara(fDx, fDy, fDz, alpha, theta, phi);
}

G4NURBS* G4Para::CreateNURBS () const
{
   // return new G4NURBSbox (fDx, fDy, fDz);
   return 0 ;
}


// ********************************  End of G4Para.cc   ********************************
