// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Cons.cc,v 1.5 1999-11-19 16:10:09 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Cons
//
// Implementation for G4Cons class
//
// History:
//
// 18.11.99 V.Grichine side = kNull initialisation in DistanceToOut(p,v,...)
// 28.04.99 V. Grichine bugs fixed in  Distance ToOut(p,v,...) and  
//          Distance ToIn(p,v)
// 09.10.98 V. Grichine modifications in Distance ToOut(p,v,...)
// 13.09.96 V. Grichine: final modifications to commit
// ~1994 P. Kent: main part of geometry functions

#include "G4Cons.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

// used by normal

enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

///////////////////////////////////////////////////////////////////////
//
// Destructor

G4Cons::~G4Cons()
{
   ;
}

//////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

G4Cons::G4Cons( const G4String& pName,
	        G4double pRmin1, G4double pRmax1,
	        G4double pRmin2, G4double pRmax2,
	        G4double pDz,
	        G4double pSPhi, G4double pDPhi) : G4CSGSolid(pName)
{
// Check z-len

    if (pDz>0)
    {
       fDz=pDz;
    }
    else
    {
       G4Exception("Error in G4Cons::G4Cons - invalid z half-length");
    }

// Check radii

    if ( pRmin1 < pRmax1 && pRmin2 < pRmax2 && pRmin1 >= 0 && pRmin2 >= 0 )
    {
       fRmin1=pRmin1; fRmax1=pRmax1;
       fRmin2=pRmin2; fRmax2=pRmax2;
    }
    else
    {
       G4Exception("Error in G4Cons::G4Cons - invalid radii");
    }

// Check angles

    if ( pDPhi >= 2.0*M_PI )
    {
       fDPhi=2*M_PI;
       fSPhi=0;
    }
    else
    {
       if ( pDPhi > 0 )
       {
	  fDPhi = pDPhi ;
       }
       else
       {
	  G4Exception("Error in G4Cons::G4Cons - invalid pDPhi");
       }
	
// Ensure pSPhi in 0-2PI or -2PI-0 range if shape crosses 0

       if (pSPhi < 0)
       {
	  fSPhi = 2.0*M_PI - fmod(fabs(pSPhi),2.0*M_PI) ;
       }
       else
       {
	  fSPhi=fmod(pSPhi,2.0*M_PI);
       }	    
       if (fSPhi+fDPhi>2.0*M_PI) fSPhi-=2.0*M_PI;
    }
}

/////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Cons::Inside(const G4ThreeVector& p) const
{
    G4double r2,rl,rh,pPhi,tolRMin,tolRMax;
    EInside in=kOutside;
    if (fabs(p.z())<=fDz-kCarTolerance/2)
	{
	    r2=p.x()*p.x()+p.y()*p.y();
	    rl=0.5*(fRmin2*(p.z()+fDz)+fRmin1*(fDz-p.z()))/fDz;
// inner radius at z of point
	    rh=0.5*(fRmax2*(p.z()+fDz)+fRmax1*(fDz-p.z()))/fDz;
// outer radius at z of point

	    if (rl)
		{
		    tolRMin=rl+kRadTolerance/2;
		}
	    else
		{
		    tolRMin=0;
		}
	    tolRMax=rh-kRadTolerance/2;
	    
	    if (r2>=tolRMin*tolRMin && r2<=tolRMax*tolRMax)
		{
		    if (fDPhi==2*M_PI||r2==0)
			{
			    in=kInside;
			}
		    else
			{
// Try inner tolerant phi boundaries (=>inside)
// if not inside, try outer tolerant phi boundaries
			    pPhi=atan2(p.y(),p.x());
			    if (pPhi<0) pPhi+=2*M_PI; // 0<=pPhi<2pi
			    if (fSPhi>=0)
				{
				    if (pPhi>=fSPhi+kAngTolerance/2 &&
					pPhi<=fSPhi+fDPhi-kAngTolerance/2)
					{
					    in=kInside;
					}
				    else if (pPhi>=fSPhi-kAngTolerance/2 &&
					     pPhi<=fSPhi+fDPhi+kAngTolerance/2)
					{
					    in=kSurface;
					}
				}
			    else
				{
				    if (pPhi<fSPhi+2*M_PI) pPhi+=2*M_PI;
				    if (pPhi>=fSPhi+2*M_PI+kAngTolerance/2 &&
					pPhi<=fSPhi+fDPhi+2*M_PI-kAngTolerance/2)
					{
					    in=kInside;
					}
				    else if (pPhi>=fSPhi+2*M_PI-kAngTolerance/2 &&
					     pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance/2)
					{
					    in=kSurface;
					}
				}
			    
			    
			}
		}
	    else
		{
// Try generous boundaries
		    tolRMin=rl-kRadTolerance/2;
		    tolRMax=rh+kRadTolerance/2;
		    if (tolRMin<0) tolRMin=0;
		    if (r2>=tolRMin*tolRMin && r2 <= tolRMax*tolRMax)
			{
			    if (fDPhi==2*M_PI||r2==0)
				{
// Continuous in phi or on z-axis
				    in=kSurface;
				}
			    else
				{
// Try outer tolerant phi boundaries only
				    pPhi=atan2(p.y(),p.x());
				    if (pPhi<0) pPhi+=2*M_PI; // 0<=pPhi<2pi
				    if (fSPhi>=0)
					{
					    if (pPhi>=fSPhi-kAngTolerance/2 &&
						pPhi<=fSPhi+fDPhi+kAngTolerance/2)
						{
						    in=kSurface;
						}
					}
				    else
					{
					    if (pPhi<fSPhi+2*M_PI) pPhi+=2*M_PI;
					    if (pPhi>=fSPhi+2*M_PI-kAngTolerance/2 &&
						pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance/2)
						{
						    in=kSurface;
						}
					}
				    

				}
			}
		}
	}
    else if (fabs(p.z())<=fDz+kCarTolerance/2)
	{
// Check within tolerant r limits
	    r2=p.x()*p.x()+p.y()*p.y();
	    rl=0.5*(fRmin2*(p.z()+fDz)+fRmin1*(fDz-p.z()))/fDz;
// inner radius at z of point
	    rh=0.5*(fRmax2*(p.z()+fDz)+fRmax1*(fDz-p.z()))/fDz;
// outer radius at z of point
	    tolRMin=rl-kRadTolerance/2;
	    if (tolRMin<0) tolRMin=0;
	    tolRMax=rh+kRadTolerance/2;
	    if (r2>=tolRMin*tolRMin && r2 <= tolRMax*tolRMax)
		{
		    if (fDPhi==2*M_PI||r2==0)
			{
// Continuous in phi or on z-axis
			    in=kSurface;
			}
		    else
			{
// Try outer tolerant phi boundaries
			    pPhi=atan2(p.y(),p.x());
			    if (pPhi<0) pPhi+=2*M_PI;		// 0<=pPhi<2pi
			    if (fSPhi>=0)
				{
				    if (pPhi>=fSPhi-kAngTolerance/2 &&
					pPhi<=fSPhi+fDPhi+kAngTolerance/2)
					{
					    in=kSurface;
					}
				}
			    else
				{
				    if (pPhi<fSPhi+2*M_PI) pPhi+=2*M_PI;
				    if (pPhi>=fSPhi+2*M_PI-kAngTolerance/2 &&
					pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance/2)
					{
					    in=kSurface;
					}
				}
			    
			}
		}
	}
    return in;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Cons::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}


///////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Cons::CalculateExtent(const EAxis pAxis,
			       const G4VoxelLimits& pVoxelLimit,
			       const G4AffineTransform& pTransform,
			       G4double& pMin,
			       G4double& pMax) const
{
    if (!pTransform.IsRotated() && fDPhi==2.0*M_PI && fRmin1==0 && fRmin2==0)
	{
// Special case handling for unrotated solid cones
// Compute z/x/y mins and maxs for bounding box respecting limits,
// with early returns if outside limits. Then switch() on pAxis,
// and compute exact x and y limit for x/y case
	    
	    G4double xoffset,xMin,xMax;
	    G4double yoffset,yMin,yMax;
	    G4double zoffset,zMin,zMax;

	    G4double diff1,diff2,maxDiff,newMin,newMax;
	    G4double xoff1,xoff2,yoff1,yoff2;
	    
	    zoffset=pTransform.NetTranslation().z();
	    zMin=zoffset-fDz;
	    zMax=zoffset+fDz;
	    if (pVoxelLimit.IsZLimited())
		{
		    if (zMin > pVoxelLimit.GetMaxZExtent()+kCarTolerance
			|| zMax < pVoxelLimit.GetMinZExtent()-kCarTolerance)
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

	    xoffset=pTransform.NetTranslation().x() ;
            G4double RMax = (fRmax2 >= fRmax1) ?  zMax : zMin  ;            	    
	    xMax = xoffset + (fRmax1 + fRmax2)*0.5 + (RMax - zoffset)*(fRmax2 - fRmax1)/(2*fDz) ;
	    xMin = 2*xoffset-xMax ;
	    if (pVoxelLimit.IsXLimited())
		{
		    if (xMin > pVoxelLimit.GetMaxXExtent()+kCarTolerance
			|| xMax < pVoxelLimit.GetMinXExtent()-kCarTolerance)
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
	    yMax = yoffset + (fRmax1 + fRmax2)*0.5 + (RMax - zoffset)*(fRmax2 - fRmax1)/(2*fDz) ;
	    yMin = 2*yoffset-yMax ;
	    RMax = yMax - yoffset ;  // is equal to max radius due to Zmax/Zmin cuttings
	    if (pVoxelLimit.IsYLimited())
		{
		    if (yMin > pVoxelLimit.GetMaxYExtent()+kCarTolerance
			|| yMax < pVoxelLimit.GetMinYExtent()-kCarTolerance)
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


// Known to cut cones
	    
	    switch (pAxis)
		{
		case kXAxis:
		    yoff1=yoffset-yMin;
		    yoff2=yMax-yoffset;
		    if (yoff1 >= 0 && yoff2 >= 0)
			{
// Y limits cross max/min x => no change
			    pMin=xMin;
			    pMax=xMax;
			}
		    else
			{
// Y limits don't cross max/min x => compute max delta x, hence new mins/maxs
			   
			    diff1=sqrt(RMax*RMax-yoff1*yoff1);
			    diff2=sqrt(RMax*RMax-yoff2*yoff2);
			    maxDiff=(diff1>diff2) ? diff1:diff2;
			    newMin=xoffset-maxDiff;
			    newMax=xoffset+maxDiff;
			    pMin=(newMin<xMin) ? xMin : newMin ;
			    pMax=(newMax>xMax) ? xMax : newMax ;
			}
	    
		    break;
		case kYAxis:
		    xoff1=xoffset-xMin;
		    xoff2=xMax-xoffset;
		    if (xoff1>=0&&xoff2>=0)
			{
// X limits cross max/min y => no change
			    pMin=yMin;
			    pMax=yMax;
			}
		    else
			{
// X limits don't cross max/min y => compute max delta y, hence new mins/maxs
			    diff1=sqrt(RMax*RMax-xoff1*xoff1);
			    diff2=sqrt(RMax*RMax-xoff2*xoff2);
			    maxDiff=(diff1>diff2) ? diff1:diff2;
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
		}

	    pMin-=kCarTolerance;
	    pMax+=kCarTolerance;

	    return true;
	    
	}
    else
	{
	    G4int i,noEntries,noBetweenSections4;
	    G4bool existsAfterClip=false;

// Calculate rotated vertex coordinates
	    G4ThreeVectorList *vertices;

	    vertices=CreateRotatedVertices(pTransform);

	    pMin=+kInfinity;
	    pMax=-kInfinity;

	    noEntries=vertices->entries();
	    noBetweenSections4=noEntries-4;
	    
	    for (i=0;i<noEntries;i+=4)
		{
		    ClipCrossSection(vertices,i,pVoxelLimit,pAxis,pMin,pMax);
		}
	    
	    for (i=0;i<noBetweenSections4;i+=4)
		{
		    ClipBetweenSections(vertices,i,pVoxelLimit,pAxis,pMin,pMax);
		}
	    
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

////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Cons::SurfaceNormal( const G4ThreeVector& p) const
{
    ENorm side;
    G4ThreeVector norm;
    G4double rho,phi;
    G4double distZ,distRMin,distRMax,distSPhi,distEPhi,distMin;
    G4double tanRMin,secRMin,pRMin,widRMin;
    G4double tanRMax,secRMax,pRMax,widRMax;

    distZ=fabs(fabs(p.z())-fDz);

    rho=sqrt(p.x()*p.x()+p.y()*p.y());

    tanRMin=(fRmin2-fRmin1)*0.5/fDz;
    secRMin=sqrt(1+tanRMin*tanRMin);
    pRMin=rho-p.z()*tanRMin;
    widRMin=fRmin2-fDz*tanRMin;
    distRMin=fabs(pRMin-widRMin)/secRMin;

    tanRMax=(fRmax2-fRmax1)*0.5/fDz;
    secRMax=sqrt(1+tanRMax*tanRMax);
    pRMax=rho-p.z()*tanRMax;
    widRMax=fRmax2-fDz*tanRMax;
    distRMax=fabs(pRMax-widRMax)/secRMax;


// First minimum
    if (distRMin<distRMax)
	{
	    if (distZ<distRMin)
		{
		    distMin=distZ;
		    side=kNZ;
		}
	    else
		{
		    distMin=distRMin;
		    side=kNRMin;
		}
	}
    else
	{
	    if (distZ<distRMax)
		{
		    distMin=distZ;
		    side=kNZ;
		}
	    else
		{
		    distMin=distRMax;
		    side=kNRMax;
		}
	}
    
    if (fDPhi<2.0*M_PI&&rho)
	{
// Protected against (0,0,z) (above)
	    phi=atan2(p.y(),p.x());
	    if (phi<0) phi+=2*M_PI;
	    if (fSPhi<0)
		{
		    distSPhi=fabs(phi-(fSPhi+2.0*M_PI))*rho;
		}
	    else
		{
		    distSPhi=fabs(phi-fSPhi)*rho;
		}

	    distEPhi=fabs(phi-fSPhi-fDPhi)*rho;

// Find new minimum
	    if (distSPhi<distEPhi)
		{
		    if (distSPhi<distMin)
			{
			    side=kNSPhi;
			}
		}
	    else
		{
		    if (distEPhi<distMin)
			{
			    side=kNEPhi;
			}
		}
	}
    
		
    switch (side)
	{
	case kNRMin:			// Inner radius
	    rho*=secRMin;
	    norm=G4ThreeVector(-p.x()/rho,-p.y()/rho,tanRMin/secRMin);
	    break;
	case kNRMax:			// Outer radius
	    rho*=secRMax;
	    norm=G4ThreeVector(p.x()/rho,p.y()/rho,-tanRMax/secRMax);
	    break;
	case kNZ:			// +/- dz
	    if (p.z()>0)
		{ norm=G4ThreeVector(0,0,1); }
	    else
		{ norm=G4ThreeVector(0,0,-1); }
	    break;
	case kNSPhi:
	    norm=G4ThreeVector(sin(fSPhi),-cos(fSPhi),0);
	    break;
	case kNEPhi:
	    norm=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
	    break;
	default:
	    G4Exception("Logic error in G4Cons::SurfaceNormal");
	    break;
	    
	} // end case

    return norm;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outside cone, compute intersection with rmax1/2
//        - if at valid phi,z return
//        - if inside outer cone, handle case when on tolerant outer cone
//          boundary and heading inwards(->0 to in)
//
// -> Compute intersection with inner cone, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points
// - z, phi intersection from Tubs

G4double G4Cons::DistanceToIn(const G4ThreeVector& p,
			      const G4ThreeVector& v) const
{
    G4double snxt=kInfinity;			// snxt = default return value

    G4bool seg;				// true if segmented in phi
    G4double hDPhi,hDPhiOT,hDPhiIT,cosHDPhiOT,cosHDPhiIT;
					// half dphi + outer tolerance
    G4double cPhi,sinCPhi,cosCPhi;	// central phi

    G4double tanRMax,secRMax,rMaxAv,rMaxOAv;	// Data for cones
    G4double tanRMin,secRMin,rMinAv,rMinIAv,rMinOAv;
    G4double rout,rin;

    G4double tolORMin,tolORMin2;	// `generous' radii squared
    G4double tolORMax2;
    G4double tolODz,tolIDz;

    G4double Dist,s,xi,yi,zi,ri,rho2,cosPsi; // Intersection point variables

    G4double t1,t2,t3,b,c,d;		// Quadratic solver variables 
    G4double nt1,nt2,nt3;
    G4double Comp;
    G4double cosSPhi,sinSPhi;		// Trig for phi start intersect
    G4double ePhi,cosEPhi,sinEPhi;	// for phi end intersect

//
// Set phi divided flag and precalcs
//
    if (fDPhi<2.0*M_PI)
	{
	    seg=true;
	    hDPhi=0.5*fDPhi;		// half delta phi
	    cPhi=fSPhi+hDPhi;;
	    hDPhiOT=hDPhi+0.5*kAngTolerance;	// outers tol' half delta phi 
	    hDPhiIT=hDPhi-0.5*kAngTolerance;
	    sinCPhi=sin(cPhi);
	    cosCPhi=cos(cPhi);
	    cosHDPhiOT=cos(hDPhiOT);
	    cosHDPhiIT=cos(hDPhiIT);
	}
    else
	{
	    seg=false;
	}

//
// Cone Precalcs
//
    tanRMin=(fRmin2-fRmin1)*0.5/fDz;
    secRMin=sqrt(1.0+tanRMin*tanRMin);
    rMinAv=(fRmin1+fRmin2)*0.5;
    if (rMinAv>kRadTolerance/2)
	{
	    rMinOAv=rMinAv-kRadTolerance/2;
	    rMinIAv=rMinAv+kRadTolerance/2;
	}
    else
	{
	    rMinOAv=0;
	    rMinIAv=0;
	}

    
    tanRMax=(fRmax2-fRmax1)*0.5/fDz;
    secRMax=sqrt(1.0+tanRMax*tanRMax);
    rMaxAv=(fRmax1+fRmax2)*0.5;

    rMaxOAv=rMaxAv+kRadTolerance/2;
   
//
// Intersection with Z surfaces
//
    tolIDz=fDz-kCarTolerance/2;
    tolODz=fDz+kCarTolerance/2;
    if (fabs(p.z())>=tolIDz)
    {
       if (p.z()*v.z()<0)		// at +Z going in -Z or visa versa
       {
	  s=(fabs(p.z())-fDz)/fabs(v.z());	// Z intersect distance

          if(s<0.0) s = 0.0 ;                  // negative dist -> zero

	  xi=p.x()+s*v.x();	// Intersection coords
	  yi=p.y()+s*v.y();
	  rho2=xi*xi+yi*yi;

          // Check validity of intersection
	  //
          // Calculate (outer) tolerant radi^2 at intersecion

	  if (v.z()>0)
	  {
	     tolORMin=fRmin1-kRadTolerance;
	     tolORMax2=(fRmax1+kRadTolerance)*(fRmax1+kRadTolerance);
	  }
	  else
	  {
	     tolORMin=fRmin2-kRadTolerance;
	     tolORMax2=(fRmax2+kRadTolerance)*(fRmax2+kRadTolerance);
	  }
	  if ( tolORMin > 0 )
	  {
			    tolORMin2=tolORMin*tolORMin;
	  }
	  else
	  {
	     tolORMin2=0;
	  }
	  if (tolORMin2 <= rho2 && rho2 <= tolORMax2)
	  {
	     if (seg&&rho2)
	     {
// Psi = angle made with central (average) phi of shape

		cosPsi=(xi*cosCPhi+yi*sinCPhi)/sqrt(rho2);

		if (cosPsi >= cosHDPhiOT)
		{
		   return s ;
		}
	     }
	     else
	     {
		return s ;
	     }
	  }
       }
       else  // On/outside extent, and heading away  -> cannot intersect
       {
	  return snxt ;	
       }
    }
    //
    // -> Can not intersect z surfaces



//
// Intersection with outer cone (possible return) and
//                   inner cone (must also check phi)
//

// Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
//
// Intersects with x^2+y^2=(a*z+b)^2
//
// where a=tanRMax or tanRMin
//       b=rMaxAv  or rMinAv
//
// (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0;
//     t1                        t2                      t3  
//
//  \--------u-------/       \-----------v----------/ \---------w--------/
//

    t1=1.0-v.z()*v.z();
    t2=p.x()*v.x()+p.y()*v.y();
    t3=p.x()*p.x()+p.y()*p.y();
    rin=tanRMin*p.z()+rMinAv;
    rout=tanRMax*p.z()+rMaxAv;

// Outer Cone Intersection
//
// Must be outside/on outer cone for valid intersection

    nt1=t1-(tanRMax*v.z())*(tanRMax*v.z());
    nt2=t2-tanRMax*v.z()*rout;
    nt3=t3-rout*rout;


    if (fabs(nt1) > kRadTolerance)  // Equation quadratic => 2 roots
    {
       if (nt3>kRadTolerance/2||rout<0)
       {
// If outside real cone (should be rho-rout>kRadTolerance/2 NOT rho^2 etc)
//                       saves a sqrt() at expense of accuracy

	  b=nt2/nt1;
	  c=nt3/nt1;
	  d=b*b-c;

          if (d>=0)
          {
			    
             if (rout < 0 && nt3 <= 0 )
             {
// Inside `shadow cone' with -ve radius -> 2nd root could be on real cone

                s = -b + sqrt(d) ;
	     }
	     else
	     {
	        if ( b <= 0 && c >= 0 ) // both >=0, try smaller
	        {
	           s = -b - sqrt(d);
	        }
                else
	        {
	           if ( c <= 0 ) // second >=0
	           {
                      s = -b + sqrt(d) ;
	           }
                   else  // both negative, travel away
	           {
                      return kInfinity ;
	           }
	        }
	     }
	     if (s>0)	// If 'forwards'. Check z intersection
	     {
	        zi=p.z()+s*v.z();

	        if (fabs(zi)<=tolODz)
	        {
// Z ok. Check phi intersection if reqd

	           if ( ! seg )
	           {
                       return s;
                   }
	           else
	           {
	              xi=p.x()+s*v.x();
	              yi=p.y()+s*v.y();
	              ri=rMaxAv+zi*tanRMax;
	              cosPsi=(xi*cosCPhi + yi*sinCPhi)/ri;

	              if (cosPsi>=cosHDPhiOT)
	              {
                         return s;
                      }
	           }
	        }             
	     }                // end if (s>0)
	  }
       }
       else
       {
// Inside outer cone
// check not inside, and heading through G4Cons (-> 0 to in)

          if (    t3 > (rin+kRadTolerance/2)*(rin+kRadTolerance/2)
               && nt2 < 0 
               && fabs(p.z()) <= tolIDz )
          {
// Inside cones, delta r -ve, inside z extent

             if (seg)
             {
        	cosPsi=(p.x()*cosCPhi+p.y()*sinCPhi)/sqrt(t3);

	        if (cosPsi>=cosHDPhiIT)
	        {
	           return 0;
	        }
	     }
	     else
	     {
	        return 0;
	     }
	  }
       }
    }
    else  //  Single root case 
    {
       if ( fabs(nt2) > kRadTolerance )
       {
          s = -0.5*nt3/nt2 ;

          if ( s < 0 )  // travel away
	  {
             return kInfinity ;
	  }
          else  // s >= 0,  If 'forwards'. Check z intersection
	  {
	     zi=p.z()+s*v.z();

	     if (fabs(zi)<=tolODz && nt2 < 0)
	     {
// Z ok. Check phi intersection if reqd

	        if ( ! seg )
	        {
                   return s;
                }
	        else
	        {
	           xi=p.x()+s*v.x();
	           yi=p.y()+s*v.y();
	           ri=rMaxAv+zi*tanRMax;
	           cosPsi=(xi*cosCPhi + yi*sinCPhi)/ri;

	           if (cosPsi>=cosHDPhiOT)
	           {
                      return s;
                   }
	        }
	     }             
	  }                
       }
       else  //    travel || cone surface from its origin
       {
          s = kInfinity ;
       }
    }


// Inner Cone Intersection
//
// o Space is divided into 3 areas:
//   1) Radius greater than real inner cone & imaginary cone & outside
//      tolerance
//   2) Radius less than inner or imaginary cone & outside tolarance
//   3) Within tolerance of real or imaginary cones
//      - Extra checks needed for 3's intersections => lots of duplicated code
    if (rMinAv)
	{
 	    nt1=t1-(tanRMin*v.z())*(tanRMin*v.z());
 	    nt2=t2-tanRMin*v.z()*rin;
 	    nt3=t3-rin*rin; 
	    if (nt1)
		{
		    if (nt3>kRadTolerance/2)
			{
// At radius greater than real & imaginary cones
// -> 2nd root, with zi check
			    b=nt2/nt1;
			    c=nt3/nt1;
			    d=b*b-c;
			    if (d>0)
				{
				    s=-b+sqrt(d);
				    if (s>0)
					{
					    zi=p.z()+s*v.z();
					    if (fabs(zi)<=tolODz)
						{
						    if (seg)
							{
 							    xi=p.x()+s*v.x();
 							    yi=p.y()+s*v.y();
							    ri=rMinAv+zi*tanRMin;
 							    cosPsi=(xi*cosCPhi+
 								    yi*sinCPhi)/ri;
 							    if (cosPsi>=cosHDPhiOT)
 								{ snxt=s; }
							}
						    else
							{
							    return s;
							}
						}
					}
				}
			    
			}
		    else  if (nt3<-kRadTolerance/2)
			{
// Within radius of inner cone (real or imaginary)
// -> Try 2nd root, with checking intersection is with real cone
// -> If check fails, try 1st root, also checking intersection is on real cone
			    b=nt2/nt1;
			    c=nt3/nt1;
			    d=b*b-c;
			    if (d>0)
				{
				    s=-b+sqrt(d);
				    zi=p.z()+s*v.z();
				    ri=rMinAv+zi*tanRMin;
				    if (ri>=0)
					{
					    if (s>0&&fabs(zi)<=tolODz)
						{
						    if (seg)
							{
							    xi=p.x()+s*v.x();
							    yi=p.y()+s*v.y();
							    cosPsi=(xi*cosCPhi+
								    yi*sinCPhi)/ri;
							    if (cosPsi>=cosHDPhiOT)
								{ snxt=s; }
							}
						    else
							{
							    return s;
							}
						}
					}
				    else
					{
					    s=-b-sqrt(d);
					    zi=p.z()+s*v.z();
					    ri=rMinAv+zi*tanRMin;
					    if (s>0&&ri>=0&&fabs(zi)<=tolODz)
						{
						    if (seg)
							{
							    xi=p.x()+s*v.x();
							    yi=p.y()+s*v.y();
							    cosPsi=(xi*cosCPhi+
								    yi*sinCPhi)/ri;
							    if (cosPsi>=cosHDPhiOT)
								{ snxt=s; }
							}
						    else
							{
							    return s;
							}
						}
					}
				}
			}
		    else
			{
// Within kRadTol/2 of inner cone (real OR imaginary)
// -> Check not travelling through (=>0 to in)
// -> if not:
//    -2nd root with validity check
//     +
			    if (fabs(p.z())<=tolODz)
				{
				    if (nt2>0)
					{
// Inside inner real cone, heading outwards, inside z range
					    if (seg)
						{
						    cosPsi=(p.x()*cosCPhi+p.y()*sinCPhi)/sqrt(t3);
						    if (cosPsi>=cosHDPhiIT)
							{
							    return 0;
							}
						}
					    else
						{
						    return 0;
						}
					}
				    else
					{
// Within z extent, but not travelling through
// -> 2nd root or kInfinity if 1st root on imaginary cone

					    b=nt2/nt1;
					    c=nt3/nt1;
					    d=b*b-c;
					    if (d>0)
						{
						    s=-b-sqrt(d);
						    zi=p.z()+s*v.z();
						    ri=rMinAv+zi*tanRMin;
						    if (ri>0)
							{
// 2nd root
							    s=-b+sqrt(d);
							    zi=p.z()+s*v.z();
							    if (s>0&&fabs(zi)<=tolODz)
								{
								    if (seg)
									{
									    xi=p.x()+s*v.x();
									    yi=p.y()+s*v.y();
									    cosPsi=(xi*cosCPhi+
										    yi*sinCPhi)/ri;
									    if (cosPsi>=cosHDPhiOT)
										{ snxt=s; }
									}
								    else
									{
									    return s;
									}
								}
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
// 2nd root
				    b=nt2/nt1;
				    c=nt3/nt1;
				    d=b*b-c;
				    if (d>0)
					{	
					    s=-b+sqrt(d);
					    zi=p.z()+s*v.z();
					    if (s>0&&fabs(zi)<=tolODz)
						{
						    if (seg)
							{
							    xi=p.x()+s*v.x();
							    yi=p.y()+s*v.y();
							    cosPsi=(xi*cosCPhi+
								    yi*sinCPhi)/ri;
							    if (cosPsi>=cosHDPhiOT)
								{ snxt=s; }
							}
						    else
							{
							    return s;
							}
						}
					}
				}
			}
		}
	}


//
// Phi segment intersection
//
// o Tolerant of points inside phi planes by up to kCarTolerance/2
//
// o NOTE: Large duplication of code between sphi & ephi checks
//         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
//            intersection check <=0 -> >=0
//         -> Should use some form of loop Construct
//
    if (seg)
	{
// First phi surface (`S'tarting phi)
	    sinSPhi=sin(fSPhi);
	    cosSPhi=cos(fSPhi);
	    Comp=v.x()*sinSPhi-v.y()*cosSPhi;
	                	// Compnent in outwards normal dirn
	    if (Comp<0)
		{
		    Dist=(p.y()*cosSPhi-p.x()*sinSPhi);
		    if (Dist<kCarTolerance/2)
			{
			    s=Dist/Comp;
			    if (s<snxt)
				{
				    if (s<0)
					{
					    s=0;
					}

				    zi=p.z()+s*v.z();
				    if (fabs(zi)<=tolODz)
					{
					    xi=p.x()+s*v.x();
					    yi=p.y()+s*v.y();
 					    rho2=xi*xi+yi*yi;
					    tolORMin2=(rMinOAv+zi*tanRMin)*(rMinOAv+zi*tanRMin);
					    tolORMax2=(rMaxOAv+zi*tanRMax)*(rMaxOAv+zi*tanRMax);

 					    if (rho2>=tolORMin2&&rho2<=tolORMax2)
 						{
// z and r intersections good - check intersecting with correct half-plane
						    if ((yi*cosCPhi-xi*sinCPhi)<=0)
							snxt=s;
 						}    
					}
				}
			}
		
		}
	    
// Second phi surface (`E'nding phi)
	    ePhi=fSPhi+fDPhi;
	    sinEPhi=sin(ePhi);
	    cosEPhi=cos(ePhi);
	    Comp=-(v.x()*sinEPhi-v.y()*cosEPhi);
				// Compnent in outwards normal dirn
	    if (Comp<0)
		{
		    Dist=-(p.y()*cosEPhi-p.x()*sinEPhi);
		    if (Dist<kCarTolerance/2)
			{
			    s=Dist/Comp;
			    if (s<snxt)
				{
				    if (s<0)
					{
					    s=0;
					}

				    zi=p.z()+s*v.z();
				    if (fabs(zi)<=tolODz)
					{
					    xi=p.x()+s*v.x();
					    yi=p.y()+s*v.y();
 					    rho2=xi*xi+yi*yi;
					    tolORMin2=(rMinOAv+zi*tanRMin)*(rMinOAv+zi*tanRMin);
					    tolORMax2=(rMaxOAv+zi*tanRMax)*(rMaxOAv+zi*tanRMax);

 					    if (rho2>=tolORMin2&&rho2<=tolORMax2)
 						{
// z and r intersections good - check intersecting with correct half-plane
						    if ((yi*cosCPhi-xi*sinCPhi)>=0)
							snxt=s;
 						}    
					}
				}
			}
		}
	}
	

	
    return snxt;
}

//////////////////////////////////////////////////////////////////////////////
// 
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Cons::DistanceToIn(const G4ThreeVector& p) const
{
    G4double safe,rho,safeR1,safeR2,safeZ;
    G4double tanRMin,secRMin,pRMin;
    G4double tanRMax,secRMax,pRMax;
    G4double phiC,cosPhiC,sinPhiC,safePhi,ePhi;
    G4double cosPsi;
    rho=sqrt(p.x()*p.x()+p.y()*p.y());
    safeZ=fabs(p.z())-fDz;

    if (fRmin1||fRmin2)
	{
	    tanRMin=(fRmin2-fRmin1)*0.5/fDz;
	    secRMin=sqrt(1.0+tanRMin*tanRMin);
	    pRMin=tanRMin*p.z()+(fRmin1+fRmin2)*0.5;
	    safeR1=(pRMin-rho)/secRMin;

	    tanRMax=(fRmax2-fRmax1)*0.5/fDz;
	    secRMax=sqrt(1.0+tanRMax*tanRMax);
	    pRMax=tanRMax*p.z()+(fRmax1+fRmax2)*0.5;
	    safeR2=(rho-pRMax)/secRMax;

	    if (safeR1>safeR2) safe=safeR1;
	    else safe=safeR2;
	}
    else
	{
	    tanRMax=(fRmax2-fRmax1)*0.5/fDz;
	    secRMax=sqrt(1.0+tanRMax*tanRMax);
	    pRMax=tanRMax*p.z()+(fRmax1+fRmax2)*0.5;
	    safe=(rho-pRMax)/secRMax;
	}
	
    if (safeZ>safe) safe=safeZ;

    if (fDPhi<2.0*M_PI&&rho)
	{
	    phiC=fSPhi+fDPhi*0.5;
	    cosPhiC=cos(phiC);
	    sinPhiC=sin(phiC);
// Psi=angle from central phi to point
	    cosPsi=(p.x()*cosPhiC+p.y()*sinPhiC)/rho;
	    if (cosPsi<cos(fDPhi*0.5))
		{
// Point lies outside phi range
		    if ((p.y()*cosPhiC-p.x()*sinPhiC)<=0)
			{
			    safePhi=fabs(p.x()*sin(fSPhi)-p.y()*cos(fSPhi));
			}
		    else
			{
			    ePhi=fSPhi+fDPhi;
			    safePhi=fabs(p.x()*sin(ePhi)-p.y()*cos(ePhi));
			}
		    if (safePhi>safe) safe=safePhi;
		}
	}
    if (safe<0) safe=0;
    return safe;
}

///////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4Cons::DistanceToOut( const G4ThreeVector& p,
			        const G4ThreeVector& v,
			        const G4bool calcNorm,
			              G4bool *validNorm,
			              G4ThreeVector *n) const
{
    ESide side = kNull,sider,sidephi;

    G4double snxt,sr,sphi,pdist;

    G4double tanRMax,secRMax,rMaxAv;	// Data for outer cone
    G4double tanRMin,secRMin,rMinAv;	// Data for inner cone

    G4double t1,t2,t3,rout,rin,nt1,nt2,nt3;
    G4double b,c,d,sr2,sr3;
    // Vars for intersection within tolerance
    ESide    sidetol;
    G4double slentol= kInfinity;
// Vars for phi intersection:
    G4double sinSPhi,cosSPhi,ePhi,sinEPhi,cosEPhi;
    G4double cPhi,sinCPhi,cosCPhi;
    G4double pDistS,compS,pDistE,compE,sphi2,xi,yi,risec,vphi;
    G4double zi,ri;
//
// Z plane intersection
//
    if (v.z()>0)
    {
       pdist=fDz-p.z();

       if (pdist > kCarTolerance*0.5)
       {
	  snxt=pdist/v.z();
	  side=kPZ;
       }
       else
       {
	  if (calcNorm)
	  {
	     *n=G4ThreeVector(0,0,1);
	     *validNorm=true;
	  }
	  return snxt=0;
       }
    }
    else if (v.z()<0)
    {
       pdist=fDz+p.z();

       if (pdist > kCarTolerance*0.5)
       {
	  snxt=-pdist/v.z();
	  side=kMZ;
       }
       else
       {
	  if (calcNorm)
	  {
	     *n=G4ThreeVector(0,0,-1);
	     *validNorm=true;
	  }
	  return snxt=0;
       }
    }
    else
    {
       snxt=kInfinity;		// Travel perpendicular to z axis
       side=kNull;
    }
//
// Radial Intersections
//
// Intersection with outer cone (possible return) and
//                   inner cone (must also check phi)
//
// Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
//
// Intersects with x^2+y^2=(a*z+b)^2
//
// where a=tanRMax or tanRMin
//       b=rMaxAv  or rMinAv
//
// (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0;
//     t1                        t2                      t3  
//
//  \--------u-------/       \-----------v----------/ \---------w--------/
//
    tanRMax=(fRmax2-fRmax1)*0.5/fDz;
    secRMax=sqrt(1.0+tanRMax*tanRMax);
    rMaxAv=(fRmax1+fRmax2)*0.5;


    t1=1.0-v.z()*v.z();			// since v normalised
    t2=p.x()*v.x()+p.y()*v.y();
    t3=p.x()*p.x()+p.y()*p.y();
    rout=tanRMax*p.z()+rMaxAv;

    nt1=t1-(tanRMax*v.z())*(tanRMax*v.z());
    nt2=t2-tanRMax*v.z()*rout;
    nt3=t3-rout*rout;
    if (nt1)
    {
//
// Equation quadratic => 2 roots : second root must be leaving
//
       b=nt2/nt1;
       c=nt3/nt1;
       d=b*b-c;

       if ( d >= 0 )
       {
// Check if on outer cone & heading outwards
// NOTE: Should use rho-rout>-kRadtolerance/2
		    
	  if (nt3 > -kRadTolerance*0.5 && nt2 >= 0 )
	  {
	     if (calcNorm)
	     {
		risec=sqrt(t3)*secRMax;
		*validNorm=true;
		*n=G4ThreeVector(p.x()/risec,p.y()/risec,-tanRMax/secRMax);
	     }
	     return snxt=0;
	  }
	  else
	 {
// // -*-* ORIG ROOT CODE    
// 			    sr=-b+sqrt(d);
// 			    sider=kRMax;
// // -*-* ORIG ROOT CODE    

// Patch 4.4.95 - root above cross-over point

	     sider=kRMax ;
	     sr=-b - sqrt(d); // was +srqrt(d), vmg 28.04.99
	     zi=p.z()+sr*v.z();
	     ri=tanRMax*zi+rMaxAv;
			    
	     if (    (ri >= 0) 
                 && (-kRadTolerance/2 <= sr)
		 && ( sr <= kRadTolerance/2)  )
	     {
// An intersection within the tolerance
//   we will Store it in case it is good -
// 
		 slentol = sr;
		 sidetol= kRMax;
	     }			      
	     if ( (ri < 0) || (sr < kRadTolerance/2) )
	     {
// Safety: if both roots -ve ensure that sr cannot `win' distance to out

		sr2=-b+sqrt(d);
		zi=p.z()+sr2*v.z();
		ri=tanRMax*zi+rMaxAv;

	        if (ri>=0&&sr2>kRadTolerance/2)
		{
		   sr=sr2;
		}
		else
		{
		   sr = kInfinity ;

		   if(    (-kRadTolerance/2 <= sr2)
		       && ( sr2 <= kRadTolerance/2)  )
		   {
// An intersection within the tolerance. Storing it in case it is good.

		      slentol = sr2;
		      sidetol= kRMax;
		   }
		}
	     }
	  }
       }
       else
       {
// No intersection with outer cone & not parallel -> already outside, no
// intersection

	  if (calcNorm)
	  {
	     risec=sqrt(t3)*secRMax;
	    *validNorm=true;
	    *n=G4ThreeVector(p.x()/risec,p.y()/risec,-tanRMax/secRMax);
	  }
	  return snxt=0;
       }
    }
    else if (nt2)
    {
//
// Linear case (only one intersection) => point outside outer cone
//
       if (calcNorm)
       {
	  risec=sqrt(t3)*secRMax;
	 *validNorm=true;
	 *n=G4ThreeVector(p.x()/risec,p.y()/risec,-tanRMax/secRMax);
       }
       return snxt=0;
    }
    else
    {
// No intersection -> parallel to outer cone => Z or inner cone intersection

       sr=kInfinity;
    }

// Check possible intersection within tolerance

    if ( slentol <= kCarTolerance/2 )
    {
//
// An intersection within the tolerance was found.  
// We must accept it only if the momentum points outwards.  
//
// G4ThreeVector ptTol;  // The point of the intersection  
// ptTol= p + slentol*v;
// ri=tanRMax*zi+rMaxAv;
//
// Calculate a normal vector,  as below

       xi=p.x()+slentol*v.x();
       yi=p.y()+slentol*v.y();
       risec=sqrt(xi*xi+yi*yi)*secRMax;
       G4ThreeVector Normal=G4ThreeVector(xi/risec,yi/risec,-tanRMax/secRMax);

       if ( Normal.dot(v) > 0 )
       {
	  // We will leave the Cone immediatelly

	  if ( calcNorm ) 
	  {
	    *n= Normal.unit();
	    *validNorm=true;
	  }
	  return snxt = 0.0;
	}
        else 
	{ 
	  // On the surface, but not heading out
	  //  so we ignore this intersection (as it is within tolerance).
	  slentol = kInfinity;
	}
    }


//
// Inner Cone intersection
//
    if (fRmin1||fRmin2)
	{
	    tanRMin=(fRmin2-fRmin1)*0.5/fDz;
	    nt1=t1-(tanRMin*v.z())*(tanRMin*v.z());
	    if (nt1)
		{
		    secRMin=sqrt(1.0+tanRMin*tanRMin);
		    rMinAv=(fRmin1+fRmin2)*0.5;
	    
		    rin=tanRMin*p.z()+rMinAv;
		    nt2=t2-tanRMin*v.z()*rin;
		    nt3=t3-rin*rin;
	    
//
// Equation quadratic => 2 roots : first root must be leaving
//
		    b=nt2/nt1;
		    c=nt3/nt1;
		    d=b*b-c;
		    if (d>=0)
			{
// NOTE: should be rho-rin<kRadTolerance/2, but using squared versions
// for efficiency
			    if (nt3<kRadTolerance*(rin+kRadTolerance*0.25)) 
				{
				    if (nt2<0)
					{
					    if (calcNorm)
						{
						    *validNorm=false;
						}
					    return snxt=0;
					}
				}
			    else
				{

// // -*-* ORIG ROOT CODE    
// 				    sr2=-b-sqrt(d);
// 				    if (sr2>0&&sr2<sr)
// 					{
// 					    sider=kRMin;
// 					    sr=sr2;
// 					}
// -*-* ORIG ROOT CODE

// // Patch 4.4.95 - zi based
// 				    sr2=-b-sqrt(d);
// 				    zi=p.z()+sr2*v.z();
// 				    if (fabs(zi)>fDz+kCarTolerance/2)
// 					{
// 					    sr2=-b+sqrt(d);
// 					    zi=p.z()+sr2*v.z();
// 					    if (fabs(zi)<fDz+kCarTolerance/2)
// 						{
// 						    if (sr2<sr&&sr2>=0)
// 							{
// 							    sr=sr2;
// 							    sider=kRMin;
// 							}
// 						}
// 					}
// 				    else
// 					{
// 					    if (sr2<sr&&sr2>=0)
// 						{
// 						    sr=sr2;
// 						    sider=kRMin;
// 						}
// 					}

// Patch 4.4.95 - root above cross-over point

				    sr2=-b-sqrt(d);
				    zi=p.z()+sr2*v.z();
				    ri=tanRMin*zi+rMinAv;
				    if( (ri>=0) 
					&& (-kRadTolerance/2 <= sr2)
					&& ( sr2 <= kRadTolerance/2)  )
				      {
				      // An intersection within the tolerance
				      //   storing it in case it is good.
					slentol = sr2;
					sidetol= kRMax;
				      }

				    if( (ri<0) 
					|| (sr2<kRadTolerance/2))
					{
					    sr3=-b+sqrt(d);

// Safety: if both roots -ve ensure that sr cannot `win' distancetoout

					    if  (sr3>kCarTolerance*0.5)
					      {
					        if(sr3<sr)
						 {
						   zi=p.z()+sr3*v.z();
						   ri=tanRMin*zi+rMinAv;
						   if (ri>=0)
						     {
						       sr=sr3;
						       sider=kRMin;
						     }
						 } 
					      }
					    else if (sr3> -kCarTolerance*0.5)
					      {
						// Intersection in tolerance
						// Store to check if it's good
						slentol = sr3;
						sidetol= kRMin;
					      }
					}
				    else if (sr2<sr&&sr2>kCarTolerance*0.5)
					{
					    sr=sr2;
					    sider=kRMin;
					}
				    else if (sr2> -kCarTolerance*0.5)
				        {
					    // Intersection in tolerance
					    // Store to check if it's good
					    slentol = sr2;
					    sidetol= kRMin;
				        }
				       
// Patch 4.4.95
				    
// // Patch 3.4.95
// 				    sr2=-b-sqrt(d);
// 				    zi=p.z()+sr2*v.z();
// 				    if (fabs(zi)>fDz+kCarTolerance/2||sr2<0)
// 					{
// 					    sr2=-b+sqrt(d);
// 					    zi=p.z()+sr2*v.z();
// 					    if (fabs(zi)>fDz+kCarTolerance/2||sr2<0)
// 						{
// 						    sr2=kInfinity;
// 						}
// 					}

//  				    if (sr2<sr)
//  					{
//  					    sider=kRMin;
//  					    sr=sr2;
// 					}
				
				    if( slentol <= kCarTolerance/2  )
				        {
					  // An intersection within the
					  //  tolerance was found. 
					  // We must accept it only if 
					  //  the momentum points outwards. 

					  G4ThreeVector Normal; 
					  
					  // Calculate a normal vector,  as below
					  xi=p.x()+slentol*v.x();
					  yi=p.y()+slentol*v.y();
					  risec=sqrt(xi*xi+yi*yi)*secRMin;
					  Normal=G4ThreeVector(xi/risec,yi/risec,-tanRMin/secRMin);
					  
					  if( Normal.dot(v) > 0 )
					    {
					      // We will leave the Cone immediatelly
					      if(calcNorm) 
						{
						  *n= Normal.unit();
						  *validNorm=true;
						}
					      return snxt = 0.0;
					    }
					  else 
					    { 
					      // On the surface, but not
					      // heading out so we ignore this
					      // intersection (as it is within
					      // tolerance). 
					      slentol = kInfinity;
					    }
					  
				        }
				}
			}
		}
//
// Linear case => point outside inner cone -> outer cone intersect
//
	}

//
// Phi Intersection
//
    
    if (fDPhi < 2.0*M_PI)
	{
	    sinSPhi=sin(fSPhi);
	    cosSPhi=cos(fSPhi);

	    ePhi=fSPhi+fDPhi;

	    sinEPhi=sin(ePhi);
	    cosEPhi=cos(ePhi);
	    cPhi=fSPhi+fDPhi*0.5;
	    sinCPhi=sin(cPhi);
	    cosCPhi=cos(cPhi);

// Check if on z axis (rho not needed later)

	    if (p.x()||p.y())
		{
// pDist -ve when inside

		    pDistS=p.x()*sinSPhi-p.y()*cosSPhi;
		    pDistE=-p.x()*sinEPhi+p.y()*cosEPhi;

// Comp -ve when in direction of outwards normal

		    compS=-sinSPhi*v.x()+cosSPhi*v.y();
		    compE=sinEPhi*v.x()-cosEPhi*v.y();

		    sidephi=kNull;

		    if (pDistS <= 0 && pDistE <= 0)
			{
// Inside both phi *full* planes
			    if (compS < 0)
				{
				    sphi=pDistS/compS;
				    xi=p.x()+sphi*v.x();
				    yi=p.y()+sphi*v.y();

// Check intersecting with correct half-plane (if not -> no intersect)
				    if ((yi*cosCPhi-xi*sinCPhi)>=0)
					{
					    sphi=kInfinity;
					}
				    else
					{
					    sidephi=kSPhi;
					    if (pDistS>-kCarTolerance/2)
						sphi=0;
					// Leave by sphi immediately
					}
				}
			    else sphi=kInfinity;
			    
			    if (compE < 0)
				{
				    sphi2=pDistE/compE;

// Only check further if < starting phi intersection

				    if (sphi2 < sphi)
					{
					    xi=p.x()+sphi2*v.x();
					    yi=p.y()+sphi2*v.y();
// Check intersecting with correct half-plane 
					    if ((yi*cosCPhi-xi*sinCPhi)>=0)
						{
// Leaving via ending phi
						    sidephi=kEPhi;

						    if (pDistE<=-kCarTolerance/2)
							{
							    sphi=sphi2;
							}
						    else 
							{
							    sphi=0;
							}
						}
					}
				}
			    
			}
		    else if (pDistS >= 0 && pDistE >= 0)
			{
// Outside both *full* phi planes
                            if (pDistS <= pDistE)
			    {
                              sidephi = kSPhi ;
			    }
                            else
			    {
                              sidephi = kEPhi ;
			    }
			    if (fDPhi>M_PI)
				{
				    if (compS < 0 && compE < 0)
                                    { 
                                        sphi=0;
				    }
				    else 
                                    {
                                        sphi=kInfinity;
				    }
				}
			    else
				{
// if towards both >=0 then once inside (after error) will remain inside

				    if (compS>=0 && compE>=0)
					{
					    sphi=kInfinity;
					}
				    else
					{
					    sphi=0;
					}
				}
			    
			}
		    else if (pDistS>0&&pDistE<0)
			{
// Outside full starting plane, inside full ending plane
			    if (fDPhi>M_PI)
				{
				    if (compE<0)
					{
					    sphi=pDistE/compE;
					    xi=p.x()+sphi*v.x();
					    yi=p.y()+sphi*v.y();
// Check intersection in correct half-plane (if not -> not leaving phi extent)
					    if ((yi*cosCPhi-xi*sinCPhi)<=0)
						{
						    sphi=kInfinity;
						}
					    else
						{
// Leaving via Ending phi
                                                    sidephi = kEPhi ;
						    if (pDistE>-kCarTolerance/2)
							sphi=0;
						}
					}
				    else
					{
					    sphi=kInfinity;
					}
				}
			    else
				{
				    if (compS>=0)
					{
					    if (compE<0)
						{
						    
						    sphi=pDistE/compE;
						    xi=p.x()+sphi*v.x();
						    yi=p.y()+sphi*v.y();
// Check intersection in correct half-plane (if not -> remain in extent)
						    if ((yi*cosCPhi-xi*sinCPhi)<=0)
							{
							    sphi=kInfinity;
							}
						    else
							{
// otherwise leaving via Ending phi
							    sidephi=kEPhi;
							}
						}
					    else sphi=kInfinity;
					}
				    else
					{
// leaving immediately by starting phi
					    sidephi=kSPhi;
					    sphi=0;
					}
				}
			}
		    else
			{
// Must be pDistS<0&&pDistE>0
// Inside full starting plane, outside full ending plane
			    if (fDPhi>M_PI)
				{
				    if (compS<0)
					{
					    sphi=pDistS/compS;
					    xi=p.x()+sphi*v.x();
					    yi=p.y()+sphi*v.y();
// Check intersection in correct half-plane (if not -> not leaving phi extent)
					    if ((yi*cosCPhi-xi*sinCPhi)>=0)
						{
						    sphi=kInfinity;
						}
					    else
						{
// Leaving via Starting phi
                                                    sidephi = kSPhi ;   
						    if (pDistS>-kCarTolerance/2)
							sphi=0;
						}
					}
				    else
					{
					    sphi=kInfinity;
					}
				}
			    else
				{
				    if (compE>=0)
					{
					    if (compS<0)
						{
						    
						    sphi=pDistS/compS;
						    xi=p.x()+sphi*v.x();
						    yi=p.y()+sphi*v.y();
// Check intersection in correct half-plane (if not -> remain in extent)
						    if ((yi*cosCPhi-xi*sinCPhi)>=0)
							    {
								sphi=kInfinity;
							    }
						    else
							{
// otherwise leaving via Starting phi
							    sidephi=kSPhi;
							}
						}
					    else
						{
						    sphi=kInfinity;
						}
					}
				    else
					{
// leaving immediately by ending
					    sidephi=kEPhi;
					    sphi=0;
					}
				}
			}
		    
		}
	    else
		{
// On z axis + travel not || to z axis -> if phi of vector direction
// within phi of shape, Step limited by rmax, else Step =0
		    vphi=atan2(v.y(),v.x());
		    if (fSPhi<vphi&&vphi<fSPhi+fDPhi)
			{
			    sphi=kInfinity;
			}
		    else
			{
                            sidephi = kSPhi ; // arbitrary 
			    sphi=0;
			}
		}
	    
// Order intersecttions
	    if (sphi<snxt)
		{
		    snxt=sphi;
		    side=sidephi;
		}
	}
// Order intersections
    if (sr<snxt)
	{
	    snxt=sr;
	    side=sider;
	}

    if (calcNorm)
	{
	    switch(side)
		{
		case kRMax:
					// Note: returned vector not normalised
					// (divide by frmax for unit vector)
		    xi=p.x()+snxt*v.x();
		    yi=p.y()+snxt*v.y();
		    risec=sqrt(xi*xi+yi*yi)*secRMax;
		    *n=G4ThreeVector(xi/risec,yi/risec,-tanRMax/secRMax);
		    *validNorm=true;
		    break;
		case kRMin:
		    *validNorm=false;	// Rmin is inconvex
		    break;
		case kSPhi:
		    if (fDPhi<=M_PI)
			{
			    *n=G4ThreeVector(sin(fSPhi),-cos(fSPhi),0);
			    *validNorm=true;
			}
		    else
			{
			    *validNorm=false;
			}
		    break;
		case kEPhi:
		    if (fDPhi<=M_PI)
			{
			*n=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
			*validNorm=true;
			}
		    else
			{
			    *validNorm=false;
			}
		    break;
		case kPZ:
		    *n=G4ThreeVector(0,0,1);
		    *validNorm=true;
		    break;
		case kMZ:
		    *n=G4ThreeVector(0,0,-1);
		    *validNorm=true;
		    break;
		default:
		    G4Exception("Invalid enum in G4Cons::DistanceToOut");
		    break;
		}
	}

    return snxt;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Cons::DistanceToOut(const G4ThreeVector& p) const
{
    G4double safe,rho,safeR1,safeR2,safeZ;
    G4double tanRMin,secRMin,pRMin;
    G4double tanRMax,secRMax,pRMax;
    G4double safePhi,phiC,cosPhiC,sinPhiC,ePhi;
    rho=sqrt(p.x()*p.x()+p.y()*p.y());

    safeZ=fDz-fabs(p.z());

    if (fRmin1||fRmin2)
	{
	    tanRMin=(fRmin2-fRmin1)*0.5/fDz;
	    secRMin=sqrt(1.0+tanRMin*tanRMin);
	    pRMin=tanRMin*p.z()+(fRmin1+fRmin2)*0.5;
	    safeR1=(rho-pRMin)/secRMin;
	}
    else
	{
	    safeR1=kInfinity;
	}

    tanRMax=(fRmax2-fRmax1)*0.5/fDz;
    secRMax=sqrt(1.0+tanRMax*tanRMax);
    pRMax=tanRMax*p.z()+(fRmax1+fRmax2)*0.5;
    safeR2=(pRMax-rho)/secRMax;

    if (safeR1<safeR2) safe=safeR1;
    else safe=safeR2;
    if (safeZ<safe) safe=safeZ;

// Check if phi divided, Calc distances closest phi plane
    if (fDPhi<2.0*M_PI)
	{
// Above/below central phi of G4Cons?
	    phiC=fSPhi+fDPhi*0.5;
	    cosPhiC=cos(phiC);
	    sinPhiC=sin(phiC);
	    if ((p.y()*cosPhiC-p.x()*sinPhiC)<=0)
		{
		    safePhi=-(p.x()*sin(fSPhi)-p.y()*cos(fSPhi));
		}
	    else
		{
		    ePhi=fSPhi+fDPhi;
		    safePhi=(p.x()*sin(ePhi)-p.y()*cos(ePhi));
		}
	    if (safePhi<safe) safe=safePhi;
	}
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
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

G4ThreeVectorList*
G4Cons::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
    G4ThreeVectorList *vertices;
    G4ThreeVector vertex0,vertex1,vertex2,vertex3;
    G4double meshAngle,meshRMax1,meshRMax2,crossAngle,cosCrossAngle,sinCrossAngle,sAngle;
    G4double rMaxX1,rMaxX2,rMaxY1,rMaxY2,rMinX1,rMinX2,rMinY1,rMinY2;
    G4int crossSection,noCrossSections;

// Compute no of cross-sections necessary to mesh cone
    
    noCrossSections=G4int (fDPhi/kMeshAngleDefault)+1;
    if (noCrossSections<kMinMeshSections)
	{
	    noCrossSections=kMinMeshSections;
	}
    else if (noCrossSections>kMaxMeshSections)
	{
	    noCrossSections=kMaxMeshSections;
	}
	
    meshAngle=fDPhi/(noCrossSections-1);

    // G4double RMax = (fRmax2 >= fRmax1) ? fRmax2 : fRmax1 ;
    meshRMax1=fRmax1/cos(meshAngle*0.5);
    meshRMax2=fRmax2/cos(meshAngle*0.5);

// If complete in phi, set start angle such that mesh will be at RMax
// on the x axis. Will give better extent calculations when not rotated.
    if (fDPhi==M_PI*2.0&&fSPhi==0)
	{
	    sAngle=-meshAngle*0.5;
	}
    else
	{
	    sAngle=fSPhi;
	}
    
    vertices=new G4ThreeVectorList(noCrossSections*4);
    if (vertices)
	{
	    for (crossSection=0;crossSection<noCrossSections;crossSection++)
		{
// Compute coordinates of cross section at section crossSection
		    crossAngle=sAngle+crossSection*meshAngle;
		    cosCrossAngle=cos(crossAngle);
		    sinCrossAngle=sin(crossAngle);

		    rMaxX1=meshRMax1*cosCrossAngle;
		    rMaxY1=meshRMax1*sinCrossAngle;
		    rMaxX2=meshRMax2*cosCrossAngle;
		    rMaxY2=meshRMax2*sinCrossAngle;
		    
		    // G4double RMin = (fRmin2 <= fRmin1) ? fRmin2 : fRmin1 ;
		    rMinX1=fRmin1*cosCrossAngle;
		    rMinY1=fRmin1*sinCrossAngle;
		    rMinX2=fRmin2*cosCrossAngle;
		    rMinY2=fRmin2*sinCrossAngle;
		    
		    vertex0=G4ThreeVector(rMinX1,rMinY1,-fDz);
		    vertex1=G4ThreeVector(rMaxX1,rMaxY1,-fDz);
		    vertex2=G4ThreeVector(rMaxX2,rMaxY2,+fDz);
		    vertex3=G4ThreeVector(rMinX2,rMinY2,+fDz);

		    vertices->insert(pTransform.TransformPoint(vertex0));
		    vertices->insert(pTransform.TransformPoint(vertex1));
		    vertices->insert(pTransform.TransformPoint(vertex2));
		    vertices->insert(pTransform.TransformPoint(vertex3));
		}
	}
    else
	{
	    G4Exception("G4Cons::CreateRotatedVertices Out of memory - Cannot alloc vertices");
	}
    return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Cons::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddThis (*this);
}

  // Define the sides of the box into which the G4Cons instance would fit.

G4VisExtent G4Cons::GetExtent() const
{
    G4double RMax = (fRmax2 >= fRmax1) ? fRmax2 : fRmax1 ;
    
    return G4VisExtent (-RMax, RMax, -RMax, RMax, -fDz, fDz);
}

G4Polyhedron* G4Cons::CreatePolyhedron () const
{
    return new G4PolyhedronCons(fRmin1, fRmax1, fRmin2, fRmax2, fDz, fSPhi,
                                fDPhi);
}

G4NURBS* G4Cons::CreateNURBS () const
{
    G4double RMax = (fRmax2 >= fRmax1) ? fRmax2 : fRmax1 ;
    
    return new G4NURBSbox (RMax, RMax, fDz);       // Box for now!!!
}

//
//
/////////////////////////////// End of G4Cons.cc file //////////////////////






