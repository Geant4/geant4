// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Trd.cc,v 1.3 1999-11-19 16:10:14 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Implementation for G4Trd class
//
// $ Id: $ 

#include "G4Trd.hh"

#include "G4VPVParameterisation.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

#include <math.h>

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths

G4Trd::G4Trd(const G4String& pName,
             G4double pdx1,  G4double pdx2,
	     G4double pdy1,  G4double pdy2,
             G4double pdz) : G4CSGSolid(pName)
{
  CheckAndSetAllParameters (pdx1, pdx2, pdy1, pdy2, pdz);
}

/////////////////////////////////////////////////////////////////////////
//
// Set and check (coplanarity) of trd parameters

void
G4Trd::CheckAndSetAllParameters (G4double pdx1,  G4double pdx2,
				 G4double pdy1,  G4double pdy2,
				 G4double pdz) 
{
  if (pdx1>0&&pdx2>0&&pdy1>0&&pdy2>0&&pdz>0)
    {
      fDx1=pdx1; fDx2=pdx2;
      fDy1=pdy1; fDy2=pdy2;
      fDz=pdz;
    }
  else
    {
      if (pdx1>=0 && pdx2>=0 && pdy1>=0 && pdy2>=0 && pdz>=0)
      {
          // G4double  Minimum_length= (1+per_thousand) * kCarTolerance/2.;
          //  FIX-ME : temporary solution for ZERO or very-small parameters.
          G4double  Minimum_length= kCarTolerance/2.;
          fDx1=max(pdx1,Minimum_length); 
          fDx2=max(pdx2,Minimum_length); 
          fDy1=max(pdy1,Minimum_length); 
          fDy2=max(pdy2,Minimum_length); 
          fDz=max(pdz,Minimum_length);
      }
      else G4Exception("Error in G4Trd::G4Trd - One or more parameters are < 0");
    }
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Trd::~G4Trd()
{}

////////////////////////////////////////////////////////////////////////////
//
//

void G4Trd::SetAllParameters (G4double pdx1, G4double pdx2, G4double pdy1, 
                              G4double pdy2, G4double pdz) 
{
  CheckAndSetAllParameters (pdx1, pdx2, pdy1, pdy2, pdz);
}


/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

 void G4Trd::ComputeDimensions(G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}


///////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Trd::CalculateExtent(const EAxis pAxis,
			      const G4VoxelLimits& pVoxelLimit,
			      const G4AffineTransform& pTransform,
			      G4double& pMin, G4double& pMax) const
{
    if (!pTransform.IsRotated())
	{
// Special case handling for unrotated solids
// Compute x/y/z mins and maxs respecting limits, with early returns
// if outside limits. Then switch() on pAxis
	    G4double xoffset,xMin,xMax;
	    G4double yoffset,yMin,yMax;
	    G4double zoffset,zMin,zMax;

	    zoffset=pTransform.NetTranslation().z();
	    zMin=zoffset-fDz;
	    zMax=zoffset+fDz;
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


	 xoffset=pTransform.NetTranslation().x();
	 if (fDx2 >= fDx1)
        { 
	 xMax =  xoffset+(fDx1+fDx2)/2+(zMax-zoffset)*(fDx2-fDx1)/(2*fDz) ;
	 xMin = 2*xoffset - xMax ;
	}
	 else
	{
	  xMax =  xoffset+(fDx1+fDx2)/2+(zMin-zoffset)*(fDx2-fDx1)/(2*fDz) ;
	  xMin =  2*xoffset - xMax ;
	}   
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


         yoffset= pTransform.NetTranslation().y() ;
	 if(fDy2 >= fDy1)
        {
	   yMax = yoffset+(fDy2+fDy1)/2+(zMax-zoffset)*(fDy2-fDy1)/(2*fDz) ;
	   yMin = 2*yoffset - yMax ;
        }
	 else
        {
	   yMax = yoffset+(fDy2+fDy1)/2+(zMin-zoffset)*(fDy2-fDy1)/(2*fDz) ;
	   yMin = 2*yoffset - yMax ;  
        }  
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
// Add 2*Tolerance to avoid precision troubles ?  
	    pMin-=kCarTolerance;
	    pMax+=kCarTolerance;

	    return true;
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
	    delete vertices;
	    return existsAfterClip;
	}

}

///////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Trd::Inside(const G4ThreeVector& p) const
{  
    EInside in=kOutside;
    double x,y,zbase1,zbase2;
    
    if (fabs(p.z())<=fDz-kCarTolerance/2)
	{
	    zbase1=p.z()+fDz;	// Dist from -ve z plane
	    zbase2=fDz-p.z();   // Dist from +ve z plane
// Check whether inside x tolerance
	    x=0.5*(fDx2*zbase1+fDx1*zbase2)/fDz - kCarTolerance/2;
	    if (fabs(p.x())<=x)
		{
		    y=0.5*((fDy2*zbase1+fDy1*zbase2))/fDz
			- kCarTolerance/2;
		    if (fabs(p.y())<=y)
			{
			    in=kInside;
			}
		    else if (fabs(p.y())<=y+kCarTolerance)
			{
			    in=kSurface;
			}
		}
	    else if (fabs(p.x())<=x+kCarTolerance)
		{
// y = y half width of shape at z of point + tolerant boundary
		    y=0.5*((fDy2*zbase1+fDy1*zbase2))/fDz
			+ kCarTolerance/2;
		    if (fabs(p.y())<=y)
			{
			    in=kSurface;
			}
		}
	}
    else if (fabs(p.z())<=fDz+kCarTolerance/2)
	{
// Only need to check outer tolerant boundaries
	    zbase1=p.z()+fDz;	// Dist from -ve z plane
	    zbase2=fDz-p.z();   // Dist from +ve z plane

// x = x half width of shape at z of point plus tolerance
	    x=0.5*(fDx2*zbase1+fDx1*zbase2)/fDz + kCarTolerance/2;
	    if (fabs(p.x())<=x)
		{
// y = y half width of shape at z of point
		    y=0.5*((fDy2*zbase1+fDy1*zbase2))/fDz
			+ kCarTolerance/2;
		    if (fabs(p.y())<=y) in=kSurface;
		}

	}	

    return in;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned

G4ThreeVector G4Trd::SurfaceNormal( const G4ThreeVector& p) const
{
   G4ThreeVector norm;
    G4double z,tanx,secx,newpx,widx;
    G4double tany,secy,newpy,widy;
    G4double distx,disty,distz,fcos;

    z=2.0*fDz;

    tanx=(fDx2-fDx1)/z;
    secx=sqrt(1.0+tanx*tanx);
    newpx=fabs(p.x())-p.z()*tanx;
    widx=fDx2-fDz*tanx;

    tany=(fDy2-fDy1)/z;
    secy=sqrt(1.0+tany*tany);
    newpy=fabs(p.y())-p.z()*tany;
    widy=fDy2-fDz*tany;

    distx=fabs(newpx-widx)/secx;	// perpendicular distance to x side
    disty=fabs(newpy-widy)/secy;	//                        to y side
    distz=fabs(fabs(p.z())-fDz);	//                        to z side

// find closest side
    if (distx<=disty)
	{ 
	    if (distx<=distz) 
		{
// Closest to X
		    fcos=1.0/secx;
//normal=(+/-cos(ang),0,-sin(ang))
		    if (p.x()>=0)
			norm=G4ThreeVector(fcos,0,-tanx*fcos);
		    else
			norm=G4ThreeVector(-fcos,0,-tanx*fcos);
		}
	    else
		{
// Closest to Z
		    if (p.z()>=0)
			norm=G4ThreeVector(0,0,1);
		    else
			norm=G4ThreeVector(0,0,-1);
		}
	}
    else
	{  
	    if (disty<=distz)
		{
// Closest to Y		    
		    fcos=1.0/secy;
		    if (p.y()>=0)
			norm=G4ThreeVector(0,fcos,-tany*fcos);
		    else
			norm=G4ThreeVector(0,-fcos,-tany*fcos);
		}
	    else 
		{
// Closest to Z
		    if (p.z()>=0)
			norm=G4ThreeVector(0,0,1);
		    else
			norm=G4ThreeVector(0,0,-1);
		}
	}
	    
    
    return norm;   
}

////////////////////////////////////////////////////////////////////////////
//
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
// NOTE:
//
// `Inside' safe - meaningful answers given if point is inside the exact
// shape.

G4double G4Trd::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v) const
{  
    G4double snxt=kInfinity;		// snxt = default return value
    G4double smin,smax;
    G4double s1,s2,tanxz,tanyz,ds1,ds2;
    G4double ss1,ss2,sn1,sn2,Dist;

// Calculate valid z intersect range
    if (v.z())
	{
// Calculate smax: must be +ve or no intersection.
	    if (v.z()>0)
		{
		    Dist=fDz-p.z();	// to plane at +dz
		    if (Dist>=kCarTolerance/2)
			{
			    smax=Dist/v.z();
			    smin=-(fDz+p.z())/v.z();
			}
		    else
			{
			    return snxt;
			}
		}
	    else
		{
// v.z <0
		    Dist=fDz+p.z();	// plane at -dz
		    if (Dist>=kCarTolerance/2)
			{
			    smax=-Dist/v.z();
			    smin=(fDz-p.z())/v.z();
			}
		    else
			{
			    return snxt;
			}
		}
	    if (smin<0) smin=0;
	}
    else
	{
// v.z=0
	    if (fabs(p.z())>fDz)
		{
		    return snxt;	// Outside & no intersect
		}
	    else
		{
		    smin=0;		// Always inside z range
		    smax=kInfinity;
		}
	}

//
// Calculate x intersection range
//

// Calc half width at p.z, and components towards planes
    tanxz=(fDx2-fDx1)*0.5/fDz;
    s1=0.5*(fDx1+fDx2)+tanxz*p.z();	// x half width at p.z
    ds1=v.x()-tanxz*v.z();		// Components of v towards faces at +-x
    ds2=v.x()+tanxz*v.z();
    ss1=s1-p.x();			// -delta x to +ve plane
					// -ve when outside
    ss2=-s1-p.x();			// -delta x to -ve plane
					// +ve when outside
    
    if (ss1<0&&ss2<=0)
	{
// In +ve coord Area
	    if (ds1<0)
		{
		    sn1=ss1/ds1;
		    if (ds2<0)
			{
			    sn2=ss2/ds2;
			}
		    else
			{
			    sn2=kInfinity;
			}
		}
	    else
		{
		    return snxt;
		}
	    
	}
    else if (ss1>=0&&ss2>0)
	{
// In -ve coord Area
	    if (ds2>0)
		{
		    sn1=ss2/ds2;
		    if (ds1>0)
			{
			    sn2=ss1/ds1;
			}
		    else
			{
			    sn2=kInfinity;
			}
		    
		}
	    else 
		{
		    return snxt;
		}
	    }
    else if (ss1>=0&&ss2<=0)
	{
// Inside Area - calculate leaving distance
// *Don't* use exact distance to side for tolerance = ss1*cos(ang xz)
//                                                  = ss1/sqrt(1.0+tanxz*tanxz)
	    sn1=0;
	    if (ds1>0)
		{
		    if (ss1>kCarTolerance/2)
			{
			    sn2=ss1/ds1; // Leave +ve side extent
			}
		    else
			{
			    return snxt; // Leave immediately by +ve
			}
		}
	    else
		{
		    sn2=kInfinity;
		}
	    
	    if (ds2<0)
		{
		    if (ss2<-kCarTolerance/2)
			{
			    Dist=ss2/ds2; // Leave -ve side extent
			    if (Dist<sn2) sn2=Dist;
			}
		    else
			{
			    return snxt;
			}
		}
	    
	}
    else if (ss1<0&&ss2>0)
	{
// Within +/- plane cross-over areas (not on boundaries ss1||ss2==0)a
	    if (ds1>=0||ds2<=0)
		{
		    return snxt;
		}
	    else
		{
// Will intersect & stay inside
		    sn1=ss1/ds1;
		    Dist=ss2/ds2;
		    if (Dist>sn1) sn1=Dist;
		    sn2=kInfinity;
		}
	}

	
// Reduce allowed range of distances as appropriate
    if (sn1>smin) smin=sn1;
    if (sn2<smax) smax=sn2;
// Check for incompatible ranges (eg z intersects between 50 ->100 and x
// only 10-40 -> no intersection)
    if (smax<smin) return snxt;


// Calculate valid y intersection range 
// (repeat of x intersection code)
    tanyz=(fDy2-fDy1)*0.5/fDz;
    s2=0.5*(fDy1+fDy2)+tanyz*p.z();	// y half width at p.z
    ds1=v.y()-tanyz*v.z();		// Components of v towards faces at +-y
    ds2=v.y()+tanyz*v.z();
    ss1=s2-p.y();			// -delta y to +ve plane
    ss2=-s2-p.y();			// -delta y to -ve plane
    
    if (ss1<0&&ss2<=0)
	{
// In +ve coord Area
	    if (ds1<0)
		{
		    sn1=ss1/ds1;
		    if (ds2<0)
			{
			    sn2=ss2/ds2;
			}
		    else
			{
			    sn2=kInfinity;
			}
		}
	    else
		{
		    return snxt;
		}
	    
	}
    else if (ss1>=0&&ss2>0)
	{
// In -ve coord Area
	    if (ds2>0)
		{
		    sn1=ss2/ds2;
		    if (ds1>0)
			{
			    sn2=ss1/ds1;
			}
		    else
			{
			    sn2=kInfinity;
			}
		    
		}
	    else 
		{
		    return snxt;
		}
	    }
    else if (ss1>=0&&ss2<=0)
	{
// Inside Area - calculate leaving distance
// *Don't* use exact distance to side for tolerance = ss1*cos(ang yz)
//                                                  = ss1/sqrt(1.0+tanyz*tanyz)
	    sn1=0;
	    if (ds1>0)
		{
		    if (ss1>kCarTolerance/2)
			{
			    sn2=ss1/ds1; // Leave +ve side extent
			}
		    else
			{
			    return snxt; // Leave immediately by +ve
			}
		}
	    else
		{
		    sn2=kInfinity;
		}
	    
	    if (ds2<0)
		{
		    if (ss2<-kCarTolerance/2)
			{
			    Dist=ss2/ds2; // Leave -ve side extent
			    if (Dist<sn2) sn2=Dist;
			}
		    else
			{
			    return snxt;
			}
		}
	    
	}
    else if (ss1<0&&ss2>0)
	{
// Within +/- plane cross-over areas (not on boundaries ss1||ss2==0)a
	    if (ds1>=0||ds2<=0)
		{
		    return snxt;
		}
	    else
		{
// Will intersect & stay inside
		    sn1=ss1/ds1;
		    Dist=ss2/ds2;
		    if (Dist>sn1) sn1=Dist;
		    sn2=kInfinity;
		}
	}
	
// Reduce allowed range of distances as appropriate
    if (sn1>smin) smin=sn1;
    if (sn2<smax) smax=sn2;

// Check for incompatible ranges (eg x intersects between 50 ->100 and y
// only 10-40 -> no intersection). Set snxt if ok
    if (smax>smin) snxt=smin;

    return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Approximate distance to shape
// Calculate perpendicular distances to z/x/y surfaces, return largest
// which is the most fast estimation of shortest distance to Trd
//  - Safe underestimate
//  - If point within exact shape, return 0 

G4double G4Trd::DistanceToIn(const G4ThreeVector& p) const
{
       G4double safe;
    G4double tanxz,distx,safx;
    G4double tanyz,disty,safy;
    G4double zbase;

    safe=fabs(p.z())-fDz;
    if (safe<0) safe=0;			// Also used to ensure x/y distances
					// POSITIVE 

    zbase=fDz+p.z();
// Find distance along x direction to closest x plane
    tanxz=(fDx2-fDx1)*0.5/fDz;
//    widx=fDx1+tanxz*(fDz+p.z()); // x width at p.z
//    distx=fabs(p.x())-widx;	   // distance to plane
    distx=fabs(p.x())-(fDx1+tanxz*zbase);
    if (distx>safe)
	{
	    safx=distx/sqrt(1.0+tanxz*tanxz); // vector Dist=Dist*cos(ang)
	    if (safx>safe) safe=safx;
	}

// Find distance along y direction to slanted wall
    tanyz=(fDy2-fDy1)*0.5/fDz;
//    widy=fDy1+tanyz*(fDz+p.z()); // y width at p.z
//    disty=fabs(p.y())-widy;	   // distance to plane
    disty=fabs(p.y())-(fDy1+tanyz*zbase);
    if (disty>safe)		
	{
	    safy=disty/sqrt(1.0+tanyz*tanyz); // distance along vector
	    if (safy>safe) safe=safy;
	}

    return safe;
}

////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance
// - z planes have std. check for tolerance
// - xz yz planes have check based on distance || to x or y axis
//   (not corrected for slope of planes)
// ?BUG? If v.z==0 are there cases when snside not set????

G4double G4Trd::DistanceToOut(const G4ThreeVector& p,
                                const G4ThreeVector& v,
			        const G4bool calcNorm,
			        G4bool *validNorm,
                                G4ThreeVector *n) const
{
    ESide side = kUndefined, snside;
    G4double snxt,pdist;
    G4double central,ss1,ss2,ds1,ds2,sn,sn2;
    G4double tanxz,cosxz,tanyz,cosyz;

    if (calcNorm) *validNorm=true; // All normals are valid

// Calculate z plane intersection
    if (v.z()>0)
	{
	    pdist=fDz-p.z();
	    if (pdist>kCarTolerance/2)
		{
		    snxt=pdist/v.z();
		    side=kPZ;
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(0,0,1);
			}
		    return snxt=0;
		}
	}
    else if (v.z()<0) 
	{
	    pdist=fDz+p.z();
	    if (pdist>kCarTolerance/2)
		{
		    snxt=-pdist/v.z();
		    side=kMZ;
		}
	    else
		{
		    if (calcNorm)
			{
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
// Calculate x intersection
//
    tanxz=(fDx2-fDx1)*0.5/fDz;
    central=0.5*(fDx1+fDx2);
// +ve plane (1)
    ss1=central+tanxz*p.z()-p.x();	// distance || x axis to plane
					// (+ve if point inside)
    ds1=v.x()-tanxz*v.z();		// component towards plane at +x
					// (-ve if +ve -> -ve direction)
// -ve plane (2)
    ss2=-tanxz*p.z()-p.x()-central;	//distance || x axis to plane
					// (-ve if point inside)
    ds2=tanxz*v.z()+v.x();		// component towards plane at -x


    if (ss1>0&&ss2<0)
	{
// Normal case - entirely inside region
	    if (ds1<=0&&ds2<0)
		{   
		    if (ss2<-kCarTolerance/2)
			{
			    sn=ss2/ds2;	// Leave by -ve side
			    snside=kMX;
			}
		    else
			{
			    sn=0; // Leave immediately by -ve side
			    snside=kMX;
			}
		}
	    else if (ds1>0&&ds2>=0)
		{
		    if (ss1>kCarTolerance/2)
			{
			    sn=ss1/ds1;	// Leave by +ve side
			    snside=kPX;
			}
		    else
			{
			    sn=0; // Leave immediately by +ve side
			    snside=kPX;
			}
		}
	    else if (ds1>0&&ds2<0)
		{
		    if (ss1>kCarTolerance/2)
			{
//				    sn=ss1/ds1;	// Leave by +ve side
			    
			    if (ss2<-kCarTolerance/2)
				{
				    sn=ss1/ds1;	// Leave by +ve side
				    sn2=ss2/ds2;
				    if (sn2<sn)
					{
					    sn=sn2;
					    snside=kMX;
					}
				    else
					{
					    snside=kPX;
					}
				}
			    else
				{
				    sn=0; // Leave immediately by -ve
				    snside=kMX;
				}
			    
			}
		    else
			{
			    sn=0; // Leave immediately by +ve side
			    snside=kPX;
			}
		    
		}
	    else
		{
// Must be || to both
		    sn=kInfinity;    // Don't leave by either side
		}
	}
    else if (ss1<=0&&ss2<0)
	{
// Outside, in +ve Area
	    if (ds1>0)
		{
		    sn=0;       // Away from shape
				        // Left by +ve side
		    snside=kPX;
		}
	    else
		{
		    if (ds2<0)
			{
// Ignore +ve plane and use -ve plane intersect
			    sn=ss2/ds2; // Leave by -ve side
			    snside=kMX;
			}
		    else
			{
// Must be || to both -> exit determined by other axes
			    sn=kInfinity; // Don't leave by either side
			}
		    
		}
	}
    else if (ss1>0&&ss2>=0)
	{
// Outside, in -ve Area
	    if (ds2<0)
		{
		    sn=0;       // away from shape
				        // Left by -ve side
		    snside=kMX;
		}
	    else
		{
		    if (ds1>0)
			{
// Ignore +ve plane and use -ve plane intersect
			    sn=ss1/ds1; // Leave by +ve side
			    snside=kPX;
			}
		    else
			{
// Must be || to both -> exit determined by other axes
			    sn=kInfinity; // Don't leave by either side
			}
		}
	}
    
// Update minimum exit distance
    if (sn<snxt)
	{
	    snxt=sn;
	    side=snside;
	}


    if (snxt>0)
	{
	    
//
// Calculate y intersection
//

	    tanyz=(fDy2-fDy1)*0.5/fDz;
	    central=0.5*(fDy1+fDy2);
// +ve plane (1)
	    ss1=central+tanyz*p.z()-p.y(); // distance || y axis to plane
					// (+ve if point inside)
	    ds1=v.y()-tanyz*v.z();	// component towards +ve plane
					// (-ve if +ve -> -ve direction)
// -ve plane (2)
	    ss2=-tanyz*p.z()-p.y()-central; // distance || y axis to plane
					// (-ve if point inside)
	    ds2=tanyz*v.z()+v.y();	// component towards -ve plane

 	    if (ss1>0&&ss2<0)
		{
// Normal case - entirely inside region
		    if (ds1<=0&&ds2<0)
			{   
			    if (ss2<-kCarTolerance/2)
				{
				    sn=ss2/ds2;	// Leave by -ve side
				    snside=kMY;
				}
			    else
				{
				    sn=0; // Leave immediately by -ve side
				    snside=kMY;
				}
			}
		    else if (ds1>0&&ds2>=0)
			{
			    if (ss1>kCarTolerance/2)
				{
				    sn=ss1/ds1;	// Leave by +ve side
				    snside=kPY;
				}
			    else
				{
				    sn=0; // Leave immediately by +ve side
				    snside=kPY;
				}
			}
		    else if (ds1>0&&ds2<0)
			{
			    if (ss1>kCarTolerance/2)
				{
//				    sn=ss1/ds1;	// Leave by +ve side
				    
				    if (ss2<-kCarTolerance/2)
					{
					    sn=ss1/ds1;	// Leave by +ve side
					    sn2=ss2/ds2;
					    if (sn2<sn)
						{
						    sn=sn2;
						    snside=kMY;
						}
					    else
						{
						    snside=kPY;
						}
					}
				    else
					{
					    sn=0; // Leave immediately by -ve
					    snside=kMY;
					}

				}
			    else
				{
				    sn=0; // Leave immediately by +ve side
				    snside=kPY;
				}
			 
			}
		    else
			{
// Must be || to both
			    sn=kInfinity;    // Don't leave by either side
			}
		}
	    else if (ss1<=0&&ss2<0)
		{
// Outside, in +ve Area
		    if (ds1>0)
			{
			    sn=0;       // Away from shape
				        // Left by +ve side
			    snside=kPY;
			}
		    else
			{
			    if (ds2<0)
				{
// Ignore +ve plane and use -ve plane intersect
				    sn=ss2/ds2; // Leave by -ve side
				    snside=kMY;
				}
			    else
				{
// Must be || to both -> exit determined by other axes
				    sn=kInfinity; // Don't leave by either side
				}

			}
		}
	    else if (ss1>0&&ss2>=0)
		{
// Outside, in -ve Area
		    if (ds2<0)
			{
			    sn=0;       // away from shape
				        // Left by -ve side
			    snside=kMY;
			}
		    else
			{
			    if (ds1>0)
				{
// Ignore +ve plane and use -ve plane intersect
				    sn=ss1/ds1; // Leave by +ve side
				    snside=kPY;
				}
			    else
				{
// Must be || to both -> exit determined by other axes
				    sn=kInfinity; // Don't leave by either side
				}
			}
		}

// Update minimum exit distance
	    if (sn<snxt)
		{
		    snxt=sn;
		    side=snside;
		}
		    
	}

    if (calcNorm)
	{
	    switch (side)
		{
		case kPX:
		    cosxz=1.0/sqrt(1.0+tanxz*tanxz);
		    *n=G4ThreeVector(cosxz,0,-tanxz*cosxz);
		    break;
		case kMX:
		    cosxz=-1.0/sqrt(1.0+tanxz*tanxz);
		    *n=G4ThreeVector(cosxz,0,tanxz*cosxz);
		    break;
		case kPY:
		    cosyz=1.0/sqrt(1.0+tanyz*tanyz);
		    *n=G4ThreeVector(0,cosyz,-tanyz*cosyz);
		    break;
		case kMY:
		    cosyz=-1.0/sqrt(1.0+tanyz*tanyz);
		    *n=G4ThreeVector(0,cosyz,tanyz*cosyz);
		    break;
		case kPZ:
		    *n=G4ThreeVector(0,0,1);
		    break;
		case kMZ:
		    *n=G4ThreeVector(0,0,-1);
		    break;
		default:
		    G4Exception("Invalid enum in G4Trd::DistanceToOut");
		    break;
		}
	}
    return snxt; 
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is point outside

G4double G4Trd::DistanceToOut(const G4ThreeVector& p) const
{
    G4double safe;
    G4double tanxz,xdist,saf1;
    G4double tanyz,ydist,saf2;
    G4double zbase;

    safe=fDz-fabs(p.z());	// z perpendicular Dist

    zbase=fDz+p.z();
// xdist = distance perpendicular to z axis to closest x plane from p
//       = (x half width of shape at p.z) - fabs(p.x)
    tanxz=(fDx2-fDx1)*0.5/fDz;
    xdist=fDx1+tanxz*zbase-fabs(p.x());
    saf1=xdist/sqrt(1.0+tanxz*tanxz); // x*cos(ang_xz) =
				// shortest (perpendicular) distance to plane
  
    tanyz=(fDy2-fDy1)*0.5/fDz;
    ydist=fDy1+tanyz*zbase-fabs(p.y());
    saf2=ydist/sqrt(1.0+tanyz*tanyz);

// Return minimum x/y/z distance
    if (safe>saf1) safe=saf1;
    if (safe>saf2) safe=saf2;

    if (safe<0) safe=0;
    return safe;	   
}

////////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility

G4ThreeVectorList*
G4Trd::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
    G4ThreeVectorList *vertices;
    vertices=new G4ThreeVectorList(8);
    if (vertices)
	{
	    G4ThreeVector vertex0(-fDx1,-fDy1,-fDz);
	    G4ThreeVector vertex1(fDx1,-fDy1,-fDz);
	    G4ThreeVector vertex2(fDx1,fDy1,-fDz);
	    G4ThreeVector vertex3(-fDx1,fDy1,-fDz);
	    G4ThreeVector vertex4(-fDx2,-fDy2,fDz);
	    G4ThreeVector vertex5(fDx2,-fDy2,fDz);
	    G4ThreeVector vertex6(fDx2,fDy2,fDz);
	    G4ThreeVector vertex7(-fDx2,fDy2,fDz);

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
	    G4Exception("G4Trd::CreateRotatedVertices Out of memory - Cannot alloc vertices");
	}
    return vertices;
}

///////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Trd::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddThis (*this);
}

G4VisExtent G4Trd::GetExtent() const
{
  G4double maxX=(fDx1>fDx2) ? fDx1 : fDx2;
  G4double maxY=(fDy1>fDy2) ? fDy1 : fDy2;
  return G4VisExtent (-maxX, maxX, -maxY, maxY, -fDz, fDz);
}

G4Polyhedron* G4Trd::CreatePolyhedron () const
{
  return new G4PolyhedronTrd2 (fDx1, fDx2, fDy1, fDy2, fDz);
}

G4NURBS* G4Trd::CreateNURBS () const
{
  //  return new G4NURBSbox (fDx, fDy, fDz);
  return 0;
}

//
//
///////////////////////////////////////////////////////////////////////////




