// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Sphere.cc,v 1.1 1999-01-07 16:07:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Sphere
//
// Implementation for G4Sphere class
//
// History:
// 28.3.94 P.Kent Old C++ code converted to tolerant geometry
// 17.9.96 V.Grichine Final modifications to commit
// 09.10.98 V. Grichine modifications in Distance ToOut(p,v,...)
// 12.11.98 V.Grichine bug fixed in DistanceToIn(p,v), theta intersections
// 25.11.98 V.Grichine bug fixed in DistanceToIn(p,v), phi intersections
// 

#include <assert.h>

#include "G4Sphere.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kSTheta,kETheta};

// used by normal

enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNSTheta,kNETheta};

// Destructor

G4Sphere::~G4Sphere()
{
   ;
}

// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4Sphere::G4Sphere(const G4String& pName,
		   G4double pRmin, G4double pRmax,
	           G4double pSPhi, G4double pDPhi,
	           G4double pSTheta, G4double pDTheta)  : G4CSGSolid(pName)
{

// Check radii
    if (pRmin<pRmax&&pRmin>=0)
	{
	   fRmin=pRmin; fRmax=pRmax;
	}
    else
	{
	    G4Exception("Error in G4Sphere::G4Sphere - invalid radii");
	}

// Check angles
    if (pDPhi>=2.0*M_PI)
	{
	   fDPhi=2*M_PI;
	}
    else if (pDPhi>0)
	{
	   fDPhi=pDPhi;
	}
    else
	{
	    G4Exception("Error in G4Sphere::G4Sphere - invalid DPhi");
	}
// Convert fSPhi to 0-2PI
    if (pSPhi<0)
	{
	   fSPhi=2.0*M_PI-fmod(fabs(pSPhi),2.0*M_PI);
	}
    else
	{
	   fSPhi=fmod(pSPhi,2.0*M_PI);
	}
// Sphere is placed such that fSPhi+fDPhi>2.0*M_PI ! fSPhi could be < 0 !!? P
    if (fSPhi+fDPhi>2.0*M_PI) fSPhi-=2.0*M_PI;

// Check theta angles
    if (pSTheta<0 || pSTheta>M_PI)
	{
	    G4Exception("Error in G4Sphere::G4Sphere - stheta outside 0-PI range");
	}
    else
	{
	   fSTheta=pSTheta;
	}

    if (pDTheta+pSTheta>=M_PI)
	{
	   fDTheta=M_PI-pSTheta;
	}
    else if (pDTheta>0)
	{
	   fDTheta=pDTheta;
	}
    else
	{
	    G4Exception("Error in G4Sphere::G4Sphere - invalid pDTheta");
	}

}

// -------------------------------------------------------------------------------------

// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
void G4Sphere::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}

// --------------------------------------------------------------------------------------

// Calculate extent under transform and specified limit
G4bool G4Sphere::CalculateExtent(const EAxis pAxis,
			      const G4VoxelLimits& pVoxelLimit,
			      const G4AffineTransform& pTransform,
			      G4double& pMin, G4double& pMax) const
{
    if ( fDPhi==2.0*M_PI && fDTheta==M_PI)  // !pTransform.IsRotated() &&
	{
// Special case handling for solid spheres-shells (rotation doesn't influence)
// Compute x/y/z mins and maxs for bounding box respecting limits,
// with early returns if outside limits. Then switch() on pAxis,
// and compute exact x and y limit for x/y case
	    
	    G4double xoffset,xMin,xMax;
	    G4double yoffset,yMin,yMax;
	    G4double zoffset,zMin,zMax;

	    G4double diff1,diff2,maxDiff,newMin,newMax;
	    G4double xoff1,xoff2,yoff1,yoff2;

	    xoffset=pTransform.NetTranslation().x();
	    xMin=xoffset-fRmax;
	    xMax=xoffset+fRmax;
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

	    yoffset=pTransform.NetTranslation().y();
	    yMin=yoffset-fRmax;
	    yMax=yoffset+fRmax;
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


	    zoffset=pTransform.NetTranslation().z();
	    zMin=zoffset-fRmax;
	    zMax=zoffset+fRmax;
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

// Known to cut sphere
	    switch (pAxis)
		{
		case kXAxis:
		    yoff1=yoffset-yMin;
		    yoff2=yMax-yoffset;
		    if (yoff1>=0&&yoff2>=0)
			{
// Y limits cross max/min x => no change
			    pMin=xMin;
			    pMax=xMax;
			}
		    else
			{
// Y limits don't cross max/min x => compute max delta x, hence new mins/maxs
			    diff1=sqrt(fRmax*fRmax-yoff1*yoff1);
			    diff2=sqrt(fRmax*fRmax-yoff2*yoff2);
			    maxDiff=(diff1>diff2) ? diff1:diff2;
			    newMin=xoffset-maxDiff;
			    newMax=xoffset+maxDiff;
			    pMin=(newMin<xMin) ? xMin : newMin;
			    pMax=(newMax>xMax) ? xMax : newMax;
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
			    diff1=sqrt(fRmax*fRmax-xoff1*xoff1);
			    diff2=sqrt(fRmax*fRmax-xoff2*xoff2);
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
    else       // Transformed cutted sphere
	{
	    G4int i,j,noEntries,noBetweenSections;
	    G4bool existsAfterClip=false;

// Calculate rotated vertex coordinates
	    G4ThreeVectorList* vertices;
            G4int  noPolygonVertices ;
	    vertices=CreateRotatedVertices(pTransform,noPolygonVertices);

	    pMin=+kInfinity;
	    pMax=-kInfinity;

	    noEntries=vertices->entries();  // noPolygonVertices*noPhiCrossSections
	    noBetweenSections=noEntries-noPolygonVertices;
	    
	         G4ThreeVectorList ThetaPolygon ;
	    for (i=0;i<noEntries;i+=noPolygonVertices)
		{
		   for(j=0;j<(noPolygonVertices/2)-1;j++)
		  {
		    ThetaPolygon.append(vertices->operator()(i+j)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+j+1)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+noPolygonVertices-2-j)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+noPolygonVertices-1-j)) ;		  
	            CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
		    ThetaPolygon.clear() ;
		  }
		}
	    for (i=0;i<noBetweenSections;i+=noPolygonVertices)
		{
		   for(j=0;j<noPolygonVertices-1;j++)
		  {
		    ThetaPolygon.append(vertices->operator()(i+j)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+j+1)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+noPolygonVertices+j+1)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+noPolygonVertices+j)) ;		  
	            CalculateClippedPolygonExtent(ThetaPolygon,pVoxelLimit,pAxis,pMin,pMax);
		    ThetaPolygon.clear() ;
		  }
		    ThetaPolygon.append(vertices->operator()(i+noPolygonVertices-1)) ;		  
		    ThetaPolygon.append(vertices->operator()(i)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+noPolygonVertices)) ;		  
		    ThetaPolygon.append(vertices->operator()(i+2*noPolygonVertices-1)) ;		  
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

// --------------------------------------------------------------------------------------------

// Return whether point inside/outside/on surface
// Split into radius, phi, theta checks
// Each check modifies `in', or returns as approprate

EInside G4Sphere::Inside(const G4ThreeVector& p) const
{
    G4double rho,rho2,rad2,tolRMin,tolRMax;
    G4double pPhi,pTheta;
    EInside in;

    rho2=p.x()*p.x()+p.y()*p.y();
    rad2=rho2+p.z()*p.z();

//
// Check radial surfaces
//  sets `in'
//
    if (fRmin)
	{
	    tolRMin=fRmin+kRadTolerance/2;
	}
    else
	{
	    tolRMin=0;
	}
    tolRMax=fRmax-kRadTolerance/2;
    
    if (rad2<=tolRMax*tolRMax&&rad2>=tolRMin*tolRMin)
	{
	    in=kInside;
	}
    else
	{
	    tolRMax=fRmax+kRadTolerance/2;
	    tolRMin=fRmin-kRadTolerance/2;
	    if (tolRMin<0)
		{
		    tolRMin=0;
		}
	    
	    if (rad2<=tolRMax*tolRMax&&rad2>=tolRMin*tolRMin)
		{
		    in=kSurface;
		}
	    else
		{
		    return in=kOutside;
		}
	}

//
// Phi boundaries   : Do not check if it has no phi boundary!
// (in!=kOutside)
//
    if ( (fDPhi<2*M_PI-kAngTolerance) && ((p.x()!=0.)||(p.y()!=0.)) )
	{
	    pPhi=atan2(p.y(),p.x());
	    if (pPhi<0) pPhi+=2*M_PI; // 0<=pPhi<2pi
	    
	    if (fSPhi>=0)
		{
		    if (in==kInside)
			{
			    if (pPhi<fSPhi+kAngTolerance/2 ||
				pPhi>fSPhi+fDPhi-kAngTolerance/2)
				{
// Not `inside' tolerant bounds
				    if (pPhi>=fSPhi-kAngTolerance/2 &&
					pPhi<=fSPhi+fDPhi+kAngTolerance/2)
					{
					    in=kSurface;
					}
				    else
					{
					    return in=kOutside;
					}
				}
			}
		    else
			{
// in==kSurface
			    if (pPhi<fSPhi-kAngTolerance/2 &&
				pPhi>fSPhi+fDPhi+kAngTolerance/2)
				{
				    return in=kOutside;
				}
			}
		}
	    else
		{
		    if (pPhi<fSPhi+2*M_PI) pPhi+=2*M_PI;
		    if (pPhi<fSPhi+2*M_PI+kAngTolerance/2 ||
			pPhi>fSPhi+fDPhi+2*M_PI-kAngTolerance/2)
			{
// Not `inside' tolerant bounds
			    if (pPhi>=fSPhi+2*M_PI-kAngTolerance/2 &&
				pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance/2)
				{
				    in=kSurface;
				}
			    else
				{
				    return in=kOutside;
				}
			}
		    
		}
	}
//
// Theta bondaries
// (in!=kOutside)
//
    if (rho2||p.z())
	{
	    rho=sqrt(rho2);
	    pTheta=atan2(rho,p.z());
	    if (in==kInside)
		{
		    if (pTheta<fSTheta+kAngTolerance/2
			|| pTheta>fSTheta+fDTheta-kAngTolerance/2)
			{
			    if (pTheta>=fSTheta-kAngTolerance/2
				&& pTheta<=fSTheta+fDTheta+kAngTolerance/2)
				{
				    in=kSurface;
				}
			    else
				{
				    in=kOutside;
				}
			}
		}
	    else
		{
		    if (pTheta<fSTheta-kAngTolerance/2
			|| pTheta>fSTheta+fDTheta+kAngTolerance/2)
			{
			    in=kOutside;
			}
		}
	}
    return in;
}

// -----------------------------------------------------------------------------------------

// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Sphere::SurfaceNormal( const G4ThreeVector& p) const
{
    ENorm side;
    G4ThreeVector norm;
    G4double rho,rho2,rad,pPhi,pTheta;
    G4double distRMin,distRMax,distSPhi,distEPhi,
	     distSTheta,distETheta,distMin;

    rho2=p.x()*p.x()+p.y()*p.y();
    rad=sqrt(rho2+p.z()*p.z());
    rho=sqrt(rho2);

//
// Distance to r shells
//

    distRMax=fabs(rad-fRmax);
    if (fRmin)
	{
	    distRMin=fabs(rad-fRmin);
	    
	    if (distRMin<distRMax)
		{
		    distMin=distRMin;
		    side=kNRMin;
		}
	    else
		{	 
		    distMin=distRMax;
		    side=kNRMax;
		}
	}
    else
	{
	    distMin=distRMax;
	    side=kNRMax;
	}

//
// Distance to phi planes
//
    if (fDPhi<2.0*M_PI&&rho)
	{
// Protected against (0,0,z) (above)
	    pPhi=atan2(p.y(),p.x());
	    if (pPhi<0) pPhi+=2*M_PI;
	    if (fSPhi<0)
		{
		    distSPhi=fabs(pPhi-(fSPhi+2.0*M_PI))*rho;
		}
	    else
		{
		    distSPhi=fabs(pPhi-fSPhi)*rho;
		}

	    distEPhi=fabs(pPhi-fSPhi-fDPhi)*rho;

// Find new minimum
	    if (distSPhi<distEPhi)
		{
		    if (distSPhi<distMin)
			{
			    distMin=distSPhi;
			    side=kNSPhi;
			}
		}
	    else
		{
		    if (distEPhi<distMin)
			{
			    distMin=distEPhi;
			    side=kNEPhi;
			}
		}
	}

//
// Distance to theta planes
//
    if (fDTheta<M_PI&&rad)
	{
	    pTheta=atan2(rho,p.z());
	    distSTheta=fabs(pTheta-fSTheta)*rad;
	    distETheta=fabs(pTheta-fSTheta-fDTheta)*rad;

// Find new minimum
	    if (distSTheta<distETheta)
		{
		    if (distSTheta<distMin)
			{
			    distMin = distSTheta ;
			    side = kNSTheta ;
			}
		}
	    else
		{
		    if (distETheta<distMin)
			{
			    distMin = distETheta ;
			    side = kNETheta ;
			}
		}
	}


    
		
    switch (side)
	{
	case kNRMin:			// Inner radius
	    norm=G4ThreeVector(-p.x()/rad,-p.y()/rad,-p.z()/rad);
	    break;
	case kNRMax:			// Outer radius
	    norm=G4ThreeVector(p.x()/rad,p.y()/rad,p.z()/rad);
	    break;
	case kNSPhi:
	    norm=G4ThreeVector(sin(fSPhi),-cos(fSPhi),0);
	    break;
	case kNEPhi:
	    norm=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
	    break;

	case kNSTheta:
	    norm=G4ThreeVector(-cos(fSTheta)*sin(fSPhi),
			       -cos(fSTheta)*cos(fSPhi),sin(fSTheta));
	    break;
	case kNETheta:
	    norm=G4ThreeVector(-cos(fSTheta+fDTheta)*cos(fSPhi+fDPhi),
			       -cos(fSTheta+fDTheta)*sin(fSPhi+fDPhi),
			       -sin(fSTheta+fSTheta));
	    break;

	default:
	    G4Exception("Logic error in G4Sphere::SurfaceNormal");
	    break;
	    
	} // end case

    return norm;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// -> If point is outside outer radius, compute intersection with rmax
//        - if no intersection return
//        - if  valid phi,theta return intersection Dist
//
// -> If shell, compute intersection with inner radius, taking largest +ve root
//        - if valid phi,theta, save intersection
//
// -> If phi segmented, compute intersection with phi half planes
//        - if valid intersection(r,theta), return smallest intersection of
//          inner shell & phi intersection
//
// -> If theta segmented, compute intersection with theta cones
//        - if valid intersection(r,phi), return smallest intersection of
//          inner shell & theta intersection
//

// NOTE:
// - `if valid' (above) implies tolerant checking of intersection points
//
// OPT:
// Move tolIO/ORmin/RMax2 precalcs to where they are needed -
// not required for most cases.
// Avoid atan2 for non theta cut G4Sphere.

G4double G4Sphere::DistanceToIn( const G4ThreeVector& p,
                                 const G4ThreeVector& v  ) const
{
    G4double snxt = kInfinity ;			// snxt = default return value

    G4double rho2, rad2, pDotV2d, pDotV3d, pTheta ;

    G4double tolIRMin2, tolORMin2, tolORMax2, tolIRMax2 ;
    G4double tolSTheta, tolETheta ;

// Intersection point

    G4double xi, yi, zi, rhoi, rhoi2, radi2, iTheta ;

// Phi intersection

    G4double sinSPhi, cosSPhi, ePhi, sinEPhi, cosEPhi , Comp ; 

// Phi flag and precalcs

    G4bool segPhi ;				
    G4double hDPhi, hDPhiOT, hDPhiIT, cPhi, sinCPhi, cosCPhi ; 
    G4double cosHDPhiOT, cosHDPhiIT ;
    G4double Dist, cosPsi ;

    G4bool segTheta ;                             // Theta flag and precals
    G4double tanSTheta, tanETheta ;
    G4double tanSTheta2, tanETheta2 ;
    G4double dist2STheta, dist2ETheta ;
    G4double t1, t2, b, c, d2, d, s = kInfinity ;

// General Precalcs

    rho2 = p.x()*p.x() + p.y()*p.y() ;
    rad2 = rho2 + p.z()*p.z() ;
    pTheta = atan2(sqrt(rho2),p.z()) ;

    pDotV2d = p.x()*v.x() + p.y()*v.y() ;
    pDotV3d = pDotV2d + p.z()*v.z() ;

// Radial Precalcs

    if (fRmin > kRadTolerance*0.5)
    {
       tolORMin2=(fRmin-kRadTolerance/2)*(fRmin-kRadTolerance/2);
    }
    else
    {
       tolORMin2 = 0 ;
    }
    tolIRMin2 = (fRmin+kRadTolerance/2)*(fRmin+kRadTolerance/2) ;
    tolORMax2 = (fRmax+kRadTolerance/2)*(fRmax+kRadTolerance/2) ;
    tolIRMax2 = (fRmax-kRadTolerance/2)*(fRmax-kRadTolerance/2) ;

// Set phi divided flag and precalcs

    if (fDPhi < 2.0*M_PI)
    {
       segPhi = true ;

       hDPhi = 0.5*fDPhi ;		// half delta phi
       cPhi = fSPhi + hDPhi ;

       hDPhiOT = hDPhi+0.5*kAngTolerance; // Outer Tolerant half delta phi 
       hDPhiIT = hDPhi-0.5*kAngTolerance;

       sinCPhi    = sin(cPhi) ;
       cosCPhi    = cos(cPhi) ;
       cosHDPhiOT = cos(hDPhiOT) ;
       cosHDPhiIT = cos(hDPhiIT) ;
    }
    else
    {
       segPhi = false ;
    }

// Theta precalcs
    
    if (fDTheta < M_PI )
    {
	segTheta  = true ;
	tolSTheta = fSTheta - kAngTolerance*0.5 ;
	tolETheta = fSTheta + fDTheta + kAngTolerance*0.5 ;
    }
    else
    {
	segTheta = false ;
    }

// Outer spherical shell intersection
// - Only if outside tolerant fRmax
// - Check for if inside and outer G4Sphere heading through solid (-> 0)
// - No intersect -> no intersection with G4Sphere
//
// Shell eqn: x^2+y^2+z^2=RSPH^2
//
// => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
//
// => (px^2+py^2+pz^2) +2s(pxvx+pyvy+pzvz)+s^2(vx^2+vy^2+vz^2)=R^2
// =>      rad2        +2s(pDotV3d)       +s^2                =R^2
//
// => s=-pDotV3d+-sqrt(pDotV3d^2-(rad2-R^2))

    c = rad2 - fRmax*fRmax ;

    if (c > kRadTolerance*fRmax)
    {

// If outside toleranct boundary of outer G4Sphere
// [should be sqrt(rad2)-fRmax > kRadTolerance/2]

	d2 = pDotV3d*pDotV3d - c ;

	if ( d2 >= 0 )
	{
	   s = -pDotV3d - sqrt(d2) ;

	   if (s >= 0 )
	   {
	      xi   = p.x() + s*v.x() ;
	      yi   = p.y() + s*v.y() ;
	      rhoi = sqrt(xi*xi + yi*yi) ;

	      if (segPhi && rhoi)    // Check phi intersection
	      {
		  cosPsi = (xi*cosCPhi + yi*sinCPhi)/rhoi ;

		  if (cosPsi >= cosHDPhiOT)
		  {
		     if (segTheta)   // Check theta intersection
		     {
			zi = p.z() + s*v.z() ;

// rhoi & zi can never both be 0 (=>intersect at origin =>fRmax=0)

			iTheta = atan2(rhoi,zi) ;

			if (iTheta >= tolSTheta && iTheta <= tolETheta)
			{
			   return snxt = s ;
			}
		     }
		     else
		     {
			return snxt=s;
		     }
		  }
	       }
	       else
	       {
		  if (segTheta)    // Check theta intersection
		  {
		     zi = p.z() + s*v.z() ;

// rhoi & zi can never both be 0 (=>intersect at origin => fRmax=0 !)

		     iTheta = atan2(rhoi,zi) ;

		     if (iTheta >= tolSTheta && iTheta <= tolETheta)
		     {
			return snxt=s;
		     }
		  }
		  else
		  {
		     return snxt = s ;
		  }
	       }			    
	    }
	 }
	 else    // No intersection with G4Sphere
	 {
	    return snxt=kInfinity;
	 }
      }
      else
      {

// Inside outer radius
// check not inside, and heading through G4Sphere (-> 0 to in)

	 if (rad2 > tolIRMin2 && pDotV3d < 0 )
	 {
	    if (segPhi)
	    {

// Use inner phi tolerant boundary -> if on tolerant
// phi boundaries, phi intersect code handles leaving/entering checks

 	       cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/sqrt(rho2) ;

	       if (cosPsi>=cosHDPhiIT)
	       { 

// inside radii, delta r -ve, inside phi

		  if (segTheta)
		  {
		     if ( pTheta >= tolSTheta + kAngTolerance && 
                          pTheta <= tolETheta - kAngTolerance     )
		     {
			return snxt=0;
		     }
		  }
		  else    // strictly inside Theta in both cases
		  {
		     return snxt=0;
		  }
	       }
	    }
	    else
	    {
	       if ( segTheta )
	       {
		  if ( pTheta >= tolSTheta + kAngTolerance && 
                       pTheta <= tolETheta - kAngTolerance     )
		  {
		     return snxt=0;
		  }
	       }
	       else   // strictly inside Theta in both cases
	       {
		  return snxt=0;
	       }
	    }
	 }

/*
//	 if (rad2 > tolIRMax2 && pDotV3d >= 0 )    // go away from convex 
//	 {                                         // outer radious surface
//	    if (segPhi)
//	    {
// Use inner phi tolerant boundary -> if on tolerant
// phi boundaries, phi intersect code handles leaving/entering checks
//
// 	        cosPsi=(p.x()*cosCPhi+p.y()*sinCPhi)/sqrt(rho2);
//
//		if (cosPsi >= cosHDPhiIT)
//		{ 
//
// inside radii, delta r -ve, inside phi
//
//		   if (segTheta)
//		   {
//		      if ( pTheta >= tolSTheta + kAngTolerance && 
//                         pTheta <= tolETheta - kAngTolerance          )
//		      {
//			 return snxt=kInfinity;
//		      }
//		   }
//		   else      // strictly inside Theta in both cases
//		   {
//		      return snxt=kInfinity;
//		   }
//		}
//	     }
//	     else
//	     {
//		if (segTheta)
//		{
//		   if ( pTheta >= tolSTheta + kAngTolerance && 
//                      pTheta <= tolETheta - kAngTolerance         )
//		   {
//		      return snxt=kInfinity;
//		   }
//		}
//		else   // strictly inside Theta in both cases
//		{
//				    return snxt=kInfinity;
//		}
//	     }
//	  }
//
  */ 

     }

// Inner spherical shell intersection
// - Always farthest root, because would have passed through outer
//   surface first.
// - Tolerant check for if travelling through solid

    if (fRmin)
    {
       c = rad2 - fRmin*fRmin ;

// Within tolerance inner radius of inner G4Sphere
// Check for immediate entry/already inside and travelling outwards

       if (c >- kRadTolerance*0.5 && pDotV3d >= 0 && rad2 < tolIRMax2 )
       {
	  if (segPhi)
	  {
// Use inner phi tolerant boundary -> if on tolerant
// phi boundaries, phi intersect code handles leaving/entering checks

	     cosPsi = (p.x()*cosCPhi+p.y()*sinCPhi)/sqrt(rho2) ;

	     if (cosPsi >= cosHDPhiIT)
	     { 

// inside radii, delta r -ve, inside phi

		 if (segTheta)
		 {
		    if ( pTheta >= tolSTheta + kAngTolerance && 
                         pTheta <= tolETheta - kAngTolerance      )
		    {
		       return snxt=0;
		    }
		 }
		 else
		 {
		    return snxt = 0 ;
		 }
	      }
	   }
	   else
	   {
	      if (segTheta)
	      {
		 if ( pTheta >= tolSTheta + kAngTolerance && 
                      pTheta <= tolETheta - kAngTolerance     )
		 {
		    return snxt = 0 ;
		 }
	      }
	      else
	      {
		 return snxt=0;
	      }
	   }
	}
	else   // Not special tolerant case
	{
	   d2 = pDotV3d*pDotV3d - c ;

	   if (d2 >= 0)
	   {
	      s = -pDotV3d + sqrt(d2) ;

	      if ( s >= kRadTolerance*0.5 )  // It was >= 0 ??
	      {
		 xi   = p.x() + s*v.x() ;
		 yi   = p.y() + s*v.y() ;
		 rhoi = sqrt(xi*xi+yi*yi) ;

		 if ( segPhi && rhoi )   // Check phi intersection
		 {
		    cosPsi = (xi*cosCPhi + yi*sinCPhi)/rhoi ;

		    if (cosPsi >= cosHDPhiOT)
		    {
		       if (segTheta)  // Check theta intersection
		       {
			  zi = p.z() + s*v.z() ;

// rhoi & zi can never both be 0 (=>intersect at origin =>fRmax=0)

			  iTheta = atan2(rhoi,zi) ;

			  if (iTheta >= tolSTheta && iTheta<=tolETheta)
			  {
			     snxt = s ;
			  }
		       }
		       else
		       {
			  snxt=s;
		       }
		    }
		 }
		 else
		 {
		    if (segTheta)   // Check theta intersection
		    {
		       zi = p.z() + s*v.z() ;

// rhoi & zi can never both be 0 (=>intersect at origin => fRmax=0 !)

		       iTheta = atan2(rhoi,zi) ;

		       if (iTheta >= tolSTheta && iTheta <= tolETheta)
		       {
			  snxt = s ;
		       }
		    }
		    else
		    {
		       snxt=s;
		    }
		 }				    
	      }

/* ******************************************************************************
//
//		// Leaving from Theta/Phi corners ?? of inner radious surface
//
//	     else if ( abs(s) < kRadTolerance && pDotV3d <= 0 )
//	     {
//		if ( segPhi )
//		{
//
// Use inner phi tolerant boundary -> if on tolerant
// phi boundaries, phi intersect code handles leaving/entering checks
//
//		   cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/sqrt(rho2) ;
//
//		   if (cosPsi<=cosHDPhiIT && cosPsi>=cosHDPhiOT)
//		   { 
//
// inside radii, delta r -ve, inside phi
//
//		      if (segTheta)
//		      {
//			 if ( pTheta >= tolSTheta + kAngTolerance && 
//                            pTheta <= tolETheta - kAngTolerance)
//			 {
//			    return snxt=kInfinity;
//			 }
//		      }
//		      else
//		      {
//			 return snxt=kInfinity;
//		      }
//		   }
//		}
//		else
//		{
//		   if ( segTheta )
//		   {
//		      if ( pTheta >= tolSTheta + kAngTolerance && 
//                         pTheta <= tolETheta - kAngTolerance     )
//		      {
//			 return snxt = kInfinity ;
//		      }
//		   }
//		   else
//		   {
//		      return snxt = kInfinity ;
//		   }
//		}
//	     }
//
   */

	  }
       }
    }

// Phi segment intersection
//
// o Tolerant of points inside phi planes by up to kCarTolerance/2
//
// o NOTE: Large duplication of code between sphi & ephi checks
//         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
//            intersection check <=0 -> >=0
//         -> Should use some form of loop Construct
//
    if ( segPhi )
    {

// First phi surface (`S'tarting phi)

       sinSPhi = sin(fSPhi) ;
       cosSPhi = cos(fSPhi) ;

// Comp = Component in outwards normal dirn

       Comp    = v.x()*sinSPhi - v.y()*cosSPhi  ;
	                	
       if ( Comp < 0 )
       {
	  Dist = p.y()*cosSPhi - p.x()*sinSPhi ;

	  if (Dist < kCarTolerance*0.5)
	  {
	     s = Dist/Comp ;

	     if (s < snxt)
	     {
		if ( s > 0 )
		{
		   xi    = p.x() + s*v.x() ;
		   yi    = p.y() + s*v.y() ;
		   zi    = p.z() + s*v.z() ;
		   rhoi2 = xi*xi + yi*yi   ;
		   radi2 = rhoi2 + zi*zi   ;
		}
	        else
		{
		   s     = 0     ;
		   xi    = p.x() ;
		   yi    = p.y() ;
		   zi    = p.z() ;
		   rhoi2 = rho2  ;
		   radi2 = rad2  ;
		}
		if ( radi2 <= tolORMax2          && 
                     radi2 >= tolORMin2          &&
                    (yi*cosCPhi-xi*sinCPhi) <= 0     )
		{

// Check theta intersection
// rhoi & zi can never both be 0 (=>intersect at origin =>fRmax=0)

		   if ( segTheta )
		   {
		      iTheta = atan2(sqrt(rhoi2),zi) ;

		      if (iTheta >= tolSTheta && iTheta <= tolETheta)
		      {
// r and theta intersections good - check intersecting with correct half-plane

			 if ((yi*cosCPhi-xi*sinCPhi) <= 0)
			 {
			    snxt = s ;
			 }
		      }
		   }
		   else
		   {
		      snxt = s ;
		   }
		}
	     }
	  }
       }

// Second phi surface (`E'nding phi)

       ePhi    = fSPhi + fDPhi ;
       sinEPhi = sin(ePhi)     ;
       cosEPhi = cos(ePhi)     ;

// Compnent in outwards normal dirn

       Comp    = -( v.x()*sinEPhi-v.y()*cosEPhi ) ;
				
       if (Comp < 0)
       {
	  Dist = -(p.y()*cosEPhi-p.x()*sinEPhi) ;

	  if ( Dist < kCarTolerance*0.5 )
	  {
	     s = Dist/Comp ;

	     if ( s < snxt )
	     {
		if (s > 0)
		{
		   xi    = p.x() + s*v.x() ;
		   yi    = p.y() + s*v.y() ;
		   zi    = p.z() + s*v.z() ;
		   rhoi2 = xi*xi + yi*yi   ;
		   radi2 = rhoi2 + zi*zi   ;
		}
		else
		{
		   s     = 0     ;
		   xi    = p.x() ;
		   yi    = p.y() ;
		   zi    = p.z() ;
		   rhoi2 = rho2  ;
		   radi2 = rad2  ;
		}
		if (  radi2 <= tolORMax2         && 
                      radi2 >= tolORMin2         &&
                     (yi*cosCPhi-xi*sinCPhi) >= 0    )
		{

// Check theta intersection
// rhoi & zi can never both be 0 (=>intersect at origin =>fRmax=0)

		   if ( segTheta )
		   {
		      iTheta = atan2(sqrt(rhoi2),zi) ;

		      if (iTheta >= tolSTheta && iTheta <= tolETheta)
		      {

// r and theta intersections good - check intersecting with correct half-plane

			 if ((yi*cosCPhi-xi*sinCPhi) >= 0)
			 {
			    snxt = s ;
			 }
		      }
		   }
		   else
		   {
		      snxt = s ;
		   }
		}
	     }
	  }
       }
    }

// Theta segment intersection

    if ( segTheta )
    {

// Intersection with theta surfaces
// Known failure cases:
// o  Inside tolerance of stheta surface, skim
//    ~parallel to cone and Hit & enter etheta surface [& visa versa]
//
//    To solve: Check 2nd root of etheta surface in addition to stheta
//
// o  start/end theta is exactly pi/2 
// Intersections with cones
//
// Cone equation: x^2+y^2=z^2tan^2(t)
//
// => (px+svx)^2+(py+svy)^2=(pz+svz)^2tan^2(t)
//
// => (px^2+py^2-pz^2tan^2(t))+2s(pxvx+pyvy-pzvztan^2(t))
//       + s^2(vx^2+vy^2-vz^2tan^2(t)) = 0
//
// => s^2(1-vz^2(1+tan^2(t))+2s(pdotv2d-pzvztan^2(t))+(rho2-pz^2tan^2(t))=0

	tanSTheta  = tan(fSTheta)         ;
	tanSTheta2 = tanSTheta*tanSTheta  ;
	tanETheta  = tan(fSTheta+fDTheta) ;
	tanETheta2 = tanETheta*tanETheta  ;
	    
	if (fSTheta)
	{
	   dist2STheta = rho2 - p.z()*p.z()*tanSTheta2 ;
	}
	else
        {
	   dist2STheta = kInfinity ;
	}
	if ( fSTheta + fDTheta < M_PI )
	{
	   dist2ETheta=rho2-p.z()*p.z()*tanETheta2;
	}
	else
	{
	   dist2ETheta=kInfinity;
	}	    
	if ( pTheta < tolSTheta)  // dist2STheta<-kRadTolerance/2 && dist2ETheta>0)
	{

// Inside (theta<stheta-tol) s theta cone
// First root of stheta cone, second if first root -ve

	   t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
	   t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;
		    
	   b  = t2/t1 ;
	   c  = dist2STheta/t1 ;
	   d2 = b*b - c ;

	   if ( d2 >= 0 )
	   {
	      d = sqrt(d2) ;
	      s = -b - d ;		// First root

	      if ( s < 0 )
	      {
		 s=-b+d;    // Second root
	      }
	      if (s >= 0 && s < snxt)
	      {
		 xi    = p.x() + s*v.x() ;
		 yi    = p.y() + s*v.y() ;
		 zi    = p.z() + s*v.z() ;
	         rhoi2 = xi*xi + yi*yi   ;
		 radi2 = rhoi2 + zi*zi   ;

		 if (radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                           zi*(fSTheta - 0.5*pi) <= 0             )
		 {
		    if ( segPhi && rhoi2 )  // Check phi intersection
		    {
		       cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

		       if (cosPsi >= cosHDPhiOT)
		       {
		          snxt = s ;
		       }
		    }
		    else
		    {
		       snxt = s ;
		    }
		 }
	      }
	   }

// Possible intersection with ETheta cone. 
// Second >= 0 root should be considered
		    
	   if ( fSTheta + fDTheta < M_PI)
	   {
	      t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
	      t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;
		    
	      b  = t2/t1 ;
	      c  = dist2ETheta/t1 ;
	      d2 = b*b - c ;

	      if (d2 >= 0)
	      {
		 d = sqrt(d2) ;
		 s = -b + d ;		// Second root

		 if (s >= 0 && s < snxt)
		 {
		    xi    = p.x() + s*v.x() ;
		    yi    = p.y() + s*v.y() ;
		    zi    = p.z() + s*v.z() ;
		    rhoi2 = xi*xi + yi*yi   ;
		    radi2 = rhoi2 + zi*zi   ;

		    if ( radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                     zi*(fSTheta + fDTheta - 0.5*pi) <= 0             )
		    {
		       if (segPhi && rhoi2)   // Check phi intersection
		       {
			  cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

			  if (cosPsi >= cosHDPhiOT)
			  {
			     snxt = s ;
			  }
		       }
		       else
		       {
			  snxt = s ;
		       }
		    }
		 }
	      }
	   }
	}  
	else if (pTheta > tolETheta) 

// dist2ETheta<-kRadTolerance/2 && dist2STheta>0)

	{
// Inside (theta>etheta+tol) e theta cone
// First root of etheta cone, second if first root `imaginary'

	   t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
	   t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;
		    
	   b  = t2/t1 ;
	   c  = dist2ETheta/t1 ;
	   d2 = b*b - c ;

	   if (d2 >= 0)
	   {
	      d = sqrt(d2) ;
	      s = -b - d ;		// First root

	      if (s < 0)
	      {
		 s = -b + d ;           // second root
	      }
	      if (s >= 0 && s < snxt)
	      {
		 xi    = p.x() + s*v.x() ;
		 yi    = p.y() + s*v.y() ;
		 zi    = p.z() + s*v.z() ;
		 rhoi2 = xi*xi + yi*yi   ;
		 radi2 = rhoi2 + zi*zi   ;

		 if (radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                 zi*(fSTheta + fDTheta - 0.5*pi) <= 0        )
		 {
		    if (segPhi && rhoi2)  // Check phi intersection
		    {
		       cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

		       if (cosPsi >= cosHDPhiOT)
		       {
			  snxt = s ;
		       }
		    }
		    else
		    {
		       snxt = s ;
		    }
		 }
	      }
	   }

// Possible intersection with STheta cone. 
// Second >= 0 root should be considered
		    
	   if ( fSTheta )
	   {
	      t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
	      t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;
		    
	      b  = t2/t1 ;
	      c  = dist2STheta/t1 ;
	      d2 = b*b - c ;

	      if (d2 >= 0)
	      {
		 d = sqrt(d2) ;
		 s = -b + d ;		// Second root

		 if (s >= 0 && s < snxt)
		 {
		    xi    = p.x() + s*v.x() ;
		    yi    = p.y() + s*v.y() ;
		    zi    = p.z() + s*v.z() ;
		    rhoi2 = xi*xi + yi*yi   ;
		    radi2 = rhoi2 + zi*zi   ;

		    if ( radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                               zi*(fSTheta - 0.5*pi) <= 0             )
		    {
		       if (segPhi && rhoi2)   // Check phi intersection
		       {
			  cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

			  if (cosPsi >= cosHDPhiOT)
			  {
			     snxt = s ;
			  }
		       }
		       else
		       {
			  snxt = s ;
		       }
		    }
		 }
	      }
	   }  
        }     
	else if ( pTheta <tolSTheta + kAngTolerance && fSTheta > kAngTolerance )  
	{

// In tolerance of stheta
// If entering through solid [r,phi] => 0 to in
// else try 2nd root

	   t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;

	   if (
              (t2>=0 && tolIRMin2<rad2 && rad2<tolIRMax2 && fSTheta<M_PI*0.5)   ||
	      (t2<0  && tolIRMin2<rad2 && rad2<tolIRMax2 && fSTheta>M_PI*0.5)   ||
              (v.z()<0 && tolIRMin2<rad2 && rad2<tolIRMax2 && fSTheta==M_PI*0.5)
	      )
	   {
	      if (segPhi && rho2)  // Check phi intersection
	      {
		 cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/sqrt(rho2) ;

		 if (cosPsi >= cosHDPhiIT)
		 {
		    return 0 ;
		 }
	      }
	      else
	      {
		 return 0 ;
	      }
	   }

// Not entering immediately/travelling through

	   t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
	   b  = t2/t1 ;
	   c  = dist2STheta/t1 ;
	   d2 = b*b - c ;

	   if (d2 >= 0)
	   {
	      d = sqrt(d2) ;
	      s = -b + d ;
	   
	      if ( s >= kCarTolerance*0.5 && s < snxt && fSTheta < M_PI*0.5 )
	      {
	         xi    = p.x() + s*v.x() ;
	         yi    = p.y() + s*v.y() ;
	         zi    = p.z() + s*v.z() ;
	         rhoi2 = xi*xi + yi*yi   ;
	         radi2 = rhoi2 + zi*zi   ;

	         if ( radi2 <= tolORMax2 && radi2 >= tolORMin2  &&
                            zi*(fSTheta - 0.5*pi) <= 0                )
	         {
	            if ( segPhi && rhoi2 )    // Check phi intersection
		    {
		       cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

		       if ( cosPsi >= cosHDPhiOT )
		       {
		          snxt = s ;
		       }
		    }
		    else
		    {
		       snxt = s ;
		    }
		 }
	      }
	   }    		    
	}   
	else if ( pTheta > tolETheta - kAngTolerance      && 
                  (fSTheta + fDTheta) < M_PI-kAngTolerance    )   
	{

// In tolerance of etheta
// If entering through solid [r,phi] => 0 to in
// else try 2nd root

	   t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;

	   if (   
    (t2<0    && (fSTheta+fDTheta) <M_PI*0.5 && tolIRMin2<rad2 && rad2<tolIRMax2)
 || (t2>=0   && (fSTheta+fDTheta) >M_PI*0.5 && tolIRMin2<rad2 && rad2<tolIRMax2)
 || (v.z()>0 && (fSTheta+fDTheta)==M_PI*0.5 && tolIRMin2<rad2 && rad2<tolIRMax2)
	      )
	   {
	      if (segPhi && rho2)   // Check phi intersection
	      {
		 cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/sqrt(rho2) ;

		 if (cosPsi >= cosHDPhiIT)
		 {
		    return 0 ;
		 }
	      }
	      else
	      {
		    return 0 ;
	      }
	   }

// Not entering immediately/travelling through

	   t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
	   b  = t2/t1 ;
	   c  = dist2ETheta/t1 ;
	   d2 = b*b - c ;

	   if (d2 >= 0)
	   {
	      d = sqrt(d2) ;
	      s = -b + d ;
	      
	      if ( s >= kCarTolerance*0.5 && 
                      s < snxt               && 
                      (fSTheta + fDTheta) > M_PI*0.5       )
	      {
		 xi    = p.x() + s*v.x() ;
		 yi    = p.y() + s*v.y() ;
		 zi    = p.z() + s*v.z() ;
		 rhoi2 = xi*xi + yi*yi   ;
		 radi2 = rhoi2 + zi*zi   ;

		 if (radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                 zi*(fSTheta + fDTheta - 0.5*pi) <= 0               )
		 {
		    if (segPhi && rhoi2)   // Check phi intersection
		    {
		       cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

		       if (cosPsi>=cosHDPhiOT)
		       {
			  snxt = s ;
		       }
		    }
		    else
		    {
		       snxt = s ;
		    }
		 }
	      }
	   }    
	}  
	else
	{

// stheta+tol<theta<etheta-tol
// For BOTH stheta & etheta check 2nd root for validity [r,phi]

	   t1 = 1 - v.z()*v.z()*(1 + tanSTheta2) ;
	   t2 = pDotV2d - p.z()*v.z()*tanSTheta2 ;

	   b  = t2/t1;
	   c  = dist2STheta/t1 ;
	   d2 = b*b - c ;

	   if (d2 >= 0)
	   {
	      d = sqrt(d2) ;
	      s = -b + d ;		// second root

	      if (s >= 0 && s < snxt)
	      {
		 xi    = p.x() + s*v.x() ;
		 yi    = p.y() + s*v.y() ;
		 zi    = p.z() + s*v.z() ;
		 rhoi2 = xi*xi + yi*yi   ;
		 radi2 = rhoi2 + zi*zi   ;

		 if (radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                     zi*(fSTheta - 0.5*pi) <= 0  )
		 {
		    if (segPhi && rhoi2)   // Check phi intersection
		    {
		       cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

		       if (cosPsi >= cosHDPhiOT)
		       {
			  snxt = s ;
		       }
		    }
		    else
		    {
		       snxt = s ;
		    }
		 }
	      }
	   }		    
	   t1 = 1 - v.z()*v.z()*(1 + tanETheta2) ;
	   t2 = pDotV2d - p.z()*v.z()*tanETheta2 ;
		    
	   b  = t2/t1 ;
	   c  = dist2ETheta/t1 ;
	   d2 = b*b - c ;

	   if (d2 >= 0)
	   {
	      d = sqrt(d2) ;
	      s = -b + d;		// second root

	      if (s >= 0 && s < snxt)
	      {
		 xi    = p.x() + s*v.x() ;
		 yi    = p.y() + s*v.y() ;
		 zi    = p.z() + s*v.z() ;
		 rhoi2 = xi*xi + yi*yi   ;
		 radi2 = rhoi2 + zi*zi   ;

		 if (radi2 <= tolORMax2 && radi2 >= tolORMin2 &&
                     zi*(fSTheta + fDTheta - 0.5*pi) <= 0           )
		 {
		    if (segPhi && rhoi2)   // Check phi intersection
		    {
		       cosPsi = (xi*cosCPhi + yi*sinCPhi)/sqrt(rhoi2) ;

		       if ( cosPsi >= cosHDPhiOT )
		       {
			  snxt=s;
		       }
		    }
		    else
		    {
		       snxt = s ;
		    }
		 }
	      }
	   }
	}  
     }

    return snxt;
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to radial planes
// - Only to phi planes if outside phi extent
// - Only to theta planes if outside theta extent
// - Return 0 if point inside

G4double G4Sphere::DistanceToIn(const G4ThreeVector& p) const
{
    G4double safe,safeRMin,safeRMax,safePhi,safeTheta;
    G4double rho2,rad,rho;
    G4double phiC,cosPhiC,sinPhiC,cosPsi,ePhi;
    G4double pTheta,dTheta1,dTheta2;
    rho2=p.x()*p.x()+p.y()*p.y();
    rad=sqrt(rho2+p.z()*p.z());
    rho=sqrt(rho2);

//
// Distance to r shells
//    
    if (fRmin)
	{
	    safeRMin=fRmin-rad;
	    safeRMax=rad-fRmax;
	    if (safeRMin>safeRMax)
		{
		    safe=safeRMin;
		}
	    else
		{
		    safe=safeRMax;
		}
	}
    else
	{
	    safe=rad-fRmax;
	}

//
// Distance to phi extent
//
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
//
// Distance to Theta extent
//    
    if ((rad!=0.0) && (fDTheta<M_PI))
	{
	    pTheta=acos(p.z()/rad);
	    if (pTheta<0) pTheta+=M_PI;
	    dTheta1=fSTheta-pTheta;
	    dTheta2=pTheta-(fSTheta+fDTheta);
	    if (dTheta1>dTheta2)
		{
		    if (dTheta1>=0)             // WHY ???????????
			{
			    safeTheta=rad*sin(dTheta1);
			    if (safe<=safeTheta)
				{
				    safe=safeTheta;
				}
			}
		}
	    else
		{
		    if (dTheta2>=0)
			{
			    safeTheta=rad*sin(dTheta2);
			    if (safe<=safeTheta)
				{
				    safe=safeTheta;
				}
			}
		    
		}
	}

    if (safe<0) safe=0;
    return safe;
}

// -----------------------------------------------------------------------------------------

// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4Sphere::DistanceToOut(const G4ThreeVector& p,
				 const G4ThreeVector& v,
			         const G4bool calcNorm,
			         G4bool *validNorm,
				 G4ThreeVector *n       ) const
{
    G4double snxt = kInfinity;     // snxt is default return value
    G4double sphi= kInfinity,stheta= kInfinity;
    ESide side=kNull,sidephi,sidetheta;  

    G4double t1,t2;
    G4double b,c,d;
                            // Vars for phi intersection:
    G4double sinSPhi,cosSPhi,ePhi,sinEPhi,cosEPhi;
    G4double cPhi,sinCPhi,cosCPhi;
    G4double pDistS,compS,pDistE,compE,sphi2,vphi;
    
    G4double rho2,rad2,pDotV2d,pDotV3d,pTheta;

    G4double tolSTheta,tolETheta;
    G4double xi,yi,zi;	    // Intersection point

    // G4double Comp; // Phi intersection

    G4bool segPhi;				// Phi flag and precalcs
    G4double hDPhi,hDPhiOT,hDPhiIT; 
    G4double cosHDPhiOT,cosHDPhiIT;

    G4bool segTheta;                             // Theta flag and precals
    G4double tanSTheta,tanETheta, rhoSecTheta;
    G4double tanSTheta2,tanETheta2;
    G4double dist2STheta,dist2ETheta;
    G4double d2,s;

//
// General Precalcs
//
    rho2=p.x()*p.x()+p.y()*p.y();
    rad2=rho2+p.z()*p.z();
    pTheta=atan2(sqrt(rho2),p.z());

    pDotV2d=p.x()*v.x()+p.y()*v.y();
    pDotV3d=pDotV2d+p.z()*v.z();


//
// Set phi divided flag and precalcs
//
    if (fDPhi<2.0*M_PI)
	{
	    segPhi=true;
	    hDPhi=0.5*fDPhi;		// half delta phi
	    cPhi=fSPhi+hDPhi;;
	    hDPhiOT=hDPhi+0.5*kAngTolerance; // Outer Tolerant half delta phi 
	    hDPhiIT=hDPhi-0.5*kAngTolerance;
	    sinCPhi=sin(cPhi);
	    cosCPhi=cos(cPhi);
	    cosHDPhiOT=cos(hDPhiOT);
	    cosHDPhiIT=cos(hDPhiIT);
	}
    else
	{
	    segPhi=false;
	}
//
// Theta precalcs
//    
    if (fDTheta<M_PI)
	{
	    segTheta=true;
	    tolSTheta=fSTheta-kAngTolerance/2;
	    tolETheta=fSTheta+fDTheta+kAngTolerance/2;
	}
    else
	{
	    segTheta=false;
	}

    
// Radial Intersections from G4Sphere::DistanceToIn
//
//
// Outer spherical shell intersection
// - Only if outside tolerant fRmax
// - Check for if inside and outer G4Sphere heading through solid (-> 0)
// - No intersect -> no intersection with G4Sphere
//
// Shell eqn: x^2+y^2+z^2=RSPH^2
//
// => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
//
// => (px^2+py^2+pz^2) +2s(pxvx+pyvy+pzvz)+s^2(vx^2+vy^2+vz^2)=R^2
// =>      rad2        +2s(pDotV3d)       +s^2                =R^2
//
// => s=-pDotV3d+-sqrt(pDotV3d^2-(rad2-R^2))
//
  G4double  Rmax_plus=  fRmax+kRadTolerance*0.5;
  G4double  Rmin_minus= (fRmin>0) ? fRmin-kRadTolerance*0.5 : 0;
  if(rad2<=Rmax_plus*Rmax_plus && rad2>=Rmin_minus*Rmin_minus)
  {
      c=rad2-fRmax*fRmax;
      if (c<kRadTolerance*fRmax) 
	{
        // Within tolerant Outer radius 
        // 
        // The test is
        //     rad  - fRmax < 0.5*kRadTolerance
        // =>  rad  < fRmax + 0.5*kRadTol
        // =>  rad2 < (fRmax + 0.5*kRadTol)^2
        // =>  rad2 < fRmax^2 + 2.*0.5*fRmax*kRadTol + 0.25*kRadTol*kRadTol
        // =>  rad2 - fRmax^2    <~    fRmax*kRadTol 

	   d2=pDotV3d*pDotV3d-c;
           if( (c>-kRadTolerance*fRmax) &&     // on tolerant surface
			       ((pDotV3d>=0)   // leaving outside from Rmax 
			       ||  (d2<0)   )) // not re-entering
	   {
	      if(calcNorm)
	      {
		 *validNorm = true ;
		 *n = G4ThreeVector(p.x()/fRmax,p.y()/fRmax,p.z()/fRmax) ;
	      }
	      return snxt = 0;
	   }
	   else 
	   {
	      snxt=-pDotV3d+sqrt(d2);    // second root since inside Rmax
	      side = kRMax ; 
	   }
	}
// Inner spherical shell intersection
// - Always first >=0 root, because would have passed from outside
//   of Rmin surface .

      if (fRmin)
	{
	    c=rad2-fRmin*fRmin;
	    if (c>-kRadTolerance*fRmin) // 2.0 * (0.5*kRadTolerance) * fRmin
		{
		   if(c<kRadTolerance*fRmin && pDotV3d<0)  // leaving from Rmin
		   {
	              if(calcNorm)
	               {
     		          *validNorm = false ;   // Rmin surface is concave
	               }
	               return snxt = 0 ;
		   }
	           else
	           {  
		    d2=pDotV3d*pDotV3d-c;
		    if (d2>=0)
		       {
			    s=-pDotV3d-sqrt(d2);
			    if (s>=0)     // Always intersect Rmin first
				{
				   snxt = s ;
				   side = kRMin ;
				}
			}
		   }
	        }
	}
  }
//
// Theta segment intersection
//
    if (segTheta)
	{
//
// Intersection with theta surfaces
//
//
// Known failure cases:
// o  Inside tolerance of stheta surface, skim
//    ~parallel to cone and Hit & enter etheta surface [& visa versa]
//
//    To solve: Check 2nd root of etheta surface in addition to stheta
//
// o  start/end theta is exactly pi/2 

//
// Intersections with cones
//
// Cone equation: x^2+y^2=z^2tan^2(t)
//
// => (px+svx)^2+(py+svy)^2=(pz+svz)^2tan^2(t)
//
// => (px^2+py^2-pz^2tan^2(t))+2s(pxvx+pyvy-pzvztan^2(t))
//       + s^2(vx^2+vy^2-vz^2tan^2(t)) = 0
//
// => s^2(1-vz^2(1+tan^2(t))+2s(pdotv2d-pzvztan^2(t))+(rho2-pz^2tan^2(t))=0
//

	    tanSTheta=tan(fSTheta);
	    tanSTheta2=tanSTheta*tanSTheta;
	    tanETheta=tan(fSTheta+fDTheta);
	    tanETheta2=tanETheta*tanETheta;
	    
	    if (fSTheta)
		{
		    dist2STheta=rho2-p.z()*p.z()*tanSTheta2;
		}
	    else
		{
		    dist2STheta=kInfinity;
		}
	    if (fSTheta+fDTheta < M_PI)
		{
		    dist2ETheta=rho2-p.z()*p.z()*tanETheta2;
		}
	    else
		{
		    dist2ETheta=kInfinity;
		}
	    
	    if (pTheta>tolSTheta && pTheta<tolETheta)   // Inside theta  
		{
// In tolerance of STheta and possible leaving out to small thetas N-
		   if(pTheta<tolSTheta+kAngTolerance  && fSTheta>kAngTolerance)  
		    {
		      t2=pDotV2d-p.z()*v.z()*tanSTheta2; // =(VdotN+)*rhoSecSTheta
		      if(fSTheta<M_PI*0.5 && t2<0)
		      {
			 if(calcNorm) *validNorm = false ;
			 return snxt = 0 ;
		      }
		     else if(fSTheta>M_PI*0.5 && t2>=0)
		      {
			 if(calcNorm)
			 {
			    rhoSecTheta = sqrt(rho2*(1+tanSTheta2)) ;
                            *validNorm = true ;
			    *n = G4ThreeVector(-p.x()/rhoSecTheta,   // N-
					       -p.y()/rhoSecTheta,
					        tanSTheta/sqrt(1+tanSTheta2)) ;
			 }
			 return snxt = 0 ;
		      }
		      else if(fSTheta==M_PI*0.5 && v.z()>0)
		      {
			 if(calcNorm)
			 {
                            *validNorm = true ;
			    *n = G4ThreeVector(0,0,1) ;
			 }
			 return snxt = 0 ;
		      }
		    }
// In tolerance of ETheta and possible leaving out to larger thetas N+
                    if(    pTheta>tolETheta-kAngTolerance 
		       && (fSTheta+fDTheta)<M_PI-kAngTolerance)  
		    {
		      t2=pDotV2d-p.z()*v.z()*tanETheta2;
		      if((fSTheta+fDTheta)>M_PI*0.5 && t2<0)
		      {
			 if(calcNorm) *validNorm = false ;
			 return snxt = 0 ;
		      }
		     else if((fSTheta+fDTheta)<M_PI*0.5 && t2>=0)
		      {
			 if(calcNorm)
			 {
			    rhoSecTheta = sqrt(rho2*(1+tanETheta2)) ;
                            *validNorm = true ;
			    *n = G4ThreeVector(p.x()/rhoSecTheta,  // N+
					       p.y()/rhoSecTheta,
					      -tanETheta/sqrt(1+tanETheta2)) ; 
			 }
			 return snxt = 0 ;
		      }
		     else if((fSTheta+fDTheta)==M_PI*0.5 && v.z()<0)
		      {
			 if(calcNorm)
			 {
                            *validNorm = true ;
			    *n = G4ThreeVector(0,0,-1) ;
			 }
			 return snxt = 0 ;
		      }
		    }
        	if(fSTheta>0)
	          {
		   
// First root of fSTheta cone, second if first root -ve
		    t1=1-v.z()*v.z()*(1+tanSTheta2);
		    t2=pDotV2d-p.z()*v.z()*tanSTheta2;
		    
		    b=t2/t1;
		    c=dist2STheta/t1;
		    d2=b*b-c;
		    if (d2>=0)
			{
			    d=sqrt(d2);
			    s=-b-d;		// First root
			    if (s<0)
				{
				    s=-b+d;    // Second root
				}
			    if (s>kRadTolerance*0.5)   // && s<sr)
				{
				   stheta = s ;
				   sidetheta = kSTheta ;
				}
			}
		  }
// Possible intersection with ETheta cone		    
	    if (fSTheta+fDTheta < M_PI)
		{
		    t1=1-v.z()*v.z()*(1+tanETheta2);
		    t2=pDotV2d-p.z()*v.z()*tanETheta2;		    
		    b=t2/t1;
		    c=dist2ETheta/t1;
		    d2=b*b-c;
		    if (d2>=0)
			{
			    d=sqrt(d2);
			    s=-b-d;	        // First root
			    if (s<0)
				{
				    s=-b+d;    // Second root
				}
			  if (s>kRadTolerance*0.5 && s<stheta)
			    {
				   stheta = s ;
				   sidetheta = kETheta ;
			    }
			}
	          }
		}  //
	}


//
// Phi Intersection
//
    
    if (fDPhi<2.0*M_PI)
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

		    if (pDistS<=0&&pDistE<=0)
			{
// Inside both phi *full* planes
			    if (compS<0)
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
			    
			    if (compE<0)
				{
				    sphi2=pDistE/compE;
// Only check further if < starting phi intersection
				    if (sphi2<sphi)
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
		    else if (pDistS>=0&&pDistE>=0)
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
				    if (compS<0&&compE<0) sphi=0;
				    else sphi=kInfinity;
				}
			    else
				{
// if towards both >=0 then once inside (after error) will remain inside
				    if (compS>=0&&compE>=0)
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

                   if ( v.x() || v.y() )
		   {
		      vphi=atan2(v.y(),v.x());

		      if ( fSPhi < vphi && vphi < fSPhi+fDPhi )
		      {
			 sphi=kInfinity;
		      }
		      else
		      {
                         sidephi = kSPhi ; // arbitrary 
			 sphi=0;
		      }
		   }
                   else  // travel along z - no phi intersaction
		   {
                      sphi = kInfinity ;
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
    if (stheta<snxt)
	{
	    snxt=stheta;
	    side=sidetheta;
	}

    if (calcNorm)    // Output switch operator
	{
	    switch(side)
		{
		case kRMax:
		    xi=p.x()+snxt*v.x();
		    yi=p.y()+snxt*v.y();
		    zi=p.z()+snxt*v.z();
		    *n=G4ThreeVector(xi/fRmax,yi/fRmax,zi/fRmax);
		    *validNorm=true;
		    break;
		case kRMin:
		    *validNorm=false;	// Rmin is concave
		    break;
		case kSPhi:
		    if (fDPhi<=M_PI)     // Normal to Phi-
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
		    if (fDPhi<=M_PI)      // Normal to Phi+
			{
			*n=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
			*validNorm=true;
			}
		    else
			{
			    *validNorm=false;
			}
		    break;
		case kSTheta:
                    if(fSTheta==M_PI*0.5)
                        {
		            *n=G4ThreeVector(0,0,1);
		            *validNorm=true;
			}
		    else if (fSTheta>M_PI)
			{
		            xi=p.x()+snxt*v.x();
		            yi=p.y()+snxt*v.y();
			    rhoSecTheta = sqrt((xi*xi+yi*yi)*(1+tanSTheta2)) ;
			    *n = G4ThreeVector(-xi/rhoSecTheta,   // N-
					       -yi/rhoSecTheta,
					        tanSTheta/sqrt(1+tanSTheta2)) ;
		            *validNorm=true;
			}
		    else
			{
			    *validNorm=false;  // Concave STheta cone
			}
		    break;
		case kETheta:
                    if((fSTheta+fDTheta)==M_PI*0.5)
                        {
		            *n=G4ThreeVector(0,0,-1);
		            *validNorm=true;
			}
		    else if ((fSTheta+fDTheta)<M_PI)
			{
		            xi=p.x()+snxt*v.x();
		            yi=p.y()+snxt*v.y();
			    rhoSecTheta = sqrt((xi*xi+yi*yi)*(1+tanETheta2)) ;
			    *n = G4ThreeVector(xi/rhoSecTheta,   // N+
					       yi/rhoSecTheta,
					      -tanSTheta/sqrt(1+tanSTheta2)) ;
		            *validNorm=true;
			}
		    else
			{
			    *validNorm=false;   // Concave ETheta cone
			}
				    break;
		default:
		    G4Exception("Invalid enum in G4Sphere::DistanceToOut");
		    break;
		}
	}
	return snxt;
}

// ----------------------------------------------------------------------------------------------

// Calcluate distance (<=actual) to closest surface of shape from inside

G4double G4Sphere::DistanceToOut(const G4ThreeVector& p) const
{
    G4double safe,safeRMin,safeRMax,safePhi,safeTheta;
    G4double rho2,rad,rho;
    G4double phiC,cosPhiC,sinPhiC,ePhi;
    G4double pTheta,dTheta1,dTheta2;
    rho2=p.x()*p.x()+p.y()*p.y();
    rad=sqrt(rho2+p.z()*p.z());
    rho=sqrt(rho2);

//
// Distance to r shells
//    
    if (fRmin)
	{
	    safeRMin=rad-fRmin;
	    safeRMax=fRmax-rad;
	    if (safeRMin<safeRMax)
		{
		    safe=safeRMin;
		}
	    else
		{
		    safe=safeRMax;
		}
	}
    else
	{
	    safe=fRmax-rad;
	}

//
// Distance to phi extent
//
    if (fDPhi<2.0*M_PI && rho)
	{
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

//
// Distance to Theta extent
//    
    if (rad)
	{
	    pTheta=acos(p.z()/rad);
	    if (pTheta<0) pTheta+=M_PI;
	    dTheta1=pTheta-fSTheta;
	    dTheta2=(fSTheta+fDTheta)-pTheta;
	    if (dTheta1<dTheta2)
		{
		    safeTheta=rad*sin(dTheta1);
		    if (safe>safeTheta)
			{
			    safe=safeTheta;
			}
		}
	    else
		{
		    safeTheta=rad*sin(dTheta2);
		    if (safe>safeTheta)
			{
			    safe=safeTheta;
			}
		}
	}

    if (safe<0) safe=0;
    return safe;
}

// -------------------------------------------------------------------------------------------

// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

G4ThreeVectorList*
G4Sphere::CreateRotatedVertices(const G4AffineTransform& pTransform,
				G4int& noPolygonVertices) const
{
    G4ThreeVectorList *vertices;
    G4ThreeVector vertex;
    G4double meshAnglePhi,meshRMax,crossAnglePhi,coscrossAnglePhi,sincrossAnglePhi,sAnglePhi;
    G4double meshTheta,crossTheta,startTheta;
    G4double rMaxX,rMaxY,rMinX,rMinY,rMinZ,rMaxZ;
    G4int crossSectionPhi,noPhiCrossSections,crossSectionTheta,noThetaSections;

// Phi cross sections
    
    noPhiCrossSections=G4int (fDPhi/kMeshAngleDefault)+1;
    
    if (noPhiCrossSections<kMinMeshSections)
	{
	    noPhiCrossSections=kMinMeshSections;
	}
    else if (noPhiCrossSections>kMaxMeshSections)
	{
	    noPhiCrossSections=kMaxMeshSections;
	}
    meshAnglePhi=fDPhi/(noPhiCrossSections-1);
    
// If complete in phi, set start angle such that mesh will be at fRMax
// on the x axis. Will give better extent calculations when not rotated.
    
    if (fDPhi==M_PI*2.0 && fSPhi==0)
	{
	    sAnglePhi = -meshAnglePhi*0.5;
	}
    else
	{
	    sAnglePhi=fSPhi;
	}    

    // Theta cross sections
    
    noThetaSections = G4int(fDTheta/kMeshAngleDefault)+1;
    
    if (noThetaSections<kMinMeshSections)
	{
	    noThetaSections=kMinMeshSections;
	}
    else if (noThetaSections>kMaxMeshSections)
	{
	    noThetaSections=kMaxMeshSections;
	}
    meshTheta=fDTheta/(noThetaSections-1);
    
// If complete in Theta, set start angle such that mesh will be at fRMax
// on the z axis. Will give better extent calculations when not rotated.
    
    if (fDTheta==M_PI && fSTheta==0)
	{
	    startTheta = -meshTheta*0.5;
	}
    else
	{
	    startTheta=fSTheta;
	}    

    meshRMax = (meshAnglePhi >= meshTheta) ?
               fRmax/cos(meshAnglePhi*0.5) : fRmax/cos(meshTheta*0.5);
    G4double* cosCrossTheta = new G4double[noThetaSections];
    G4double* sinCrossTheta = new G4double[noThetaSections];    
    vertices=new G4ThreeVectorList(noPhiCrossSections*(noThetaSections*2));
    if (vertices && cosCrossTheta && sinCrossTheta)
	{
    for (crossSectionPhi=0;crossSectionPhi<noPhiCrossSections;crossSectionPhi++)
          {
		  crossAnglePhi=sAnglePhi+crossSectionPhi*meshAnglePhi;
		  coscrossAnglePhi=cos(crossAnglePhi);
		  sincrossAnglePhi=sin(crossAnglePhi);
    for (crossSectionTheta=0;crossSectionTheta<noThetaSections;crossSectionTheta++)
		{
// Compute coordinates of cross section at section crossSectionPhi
		   
		  crossTheta=startTheta+crossSectionTheta*meshTheta;
		  cosCrossTheta[crossSectionTheta]=cos(crossTheta);
		  sinCrossTheta[crossSectionTheta]=sin(crossTheta);

		  rMinX=fRmin*sinCrossTheta[crossSectionTheta]*coscrossAnglePhi;
		  rMinY=fRmin*sinCrossTheta[crossSectionTheta]*sincrossAnglePhi;
		  rMinZ=fRmin*cosCrossTheta[crossSectionTheta];
		    
		    vertex=G4ThreeVector(rMinX,rMinY,rMinZ);
		    vertices->insert(pTransform.TransformPoint(vertex));
		    
		}    // Theta forward 
    
    for (crossSectionTheta=noThetaSections-1;crossSectionTheta>=0;crossSectionTheta--)
		{
		  rMaxX=meshRMax*sinCrossTheta[crossSectionTheta]*coscrossAnglePhi;
		  rMaxY=meshRMax*sinCrossTheta[crossSectionTheta]*sincrossAnglePhi;
		  rMaxZ=meshRMax*cosCrossTheta[crossSectionTheta];
		    
		    vertex=G4ThreeVector(rMaxX,rMaxY,rMaxZ);
		    vertices->insert(pTransform.TransformPoint(vertex));

		}   // Theta back 
	  }       // Phi
                                   noPolygonVertices = noThetaSections*2 ;
	}
    else
	{
	    G4Exception("G4Sphere::CreateRotatedVertices Out of memory - Cannot alloc vertices");
	}

    delete[] cosCrossTheta;
    delete[] sinCrossTheta;

    return vertices;
}

// ---------------------------------------------------------------------------------------

void G4Sphere::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddThis (*this);
}

G4VisExtent G4Sphere::GetExtent() const
{
  // Define the sides of the box into which the G4Sphere instance would fit.
  return G4VisExtent (-fRmax, fRmax, -fRmax, fRmax, -fRmax, fRmax);
}

G4Polyhedron* G4Sphere::CreatePolyhedron () const
{
    return new G4PolyhedronSphere (fRmin, fRmax,
				   fSPhi, fDPhi,
				   fSTheta, fDTheta);

}

G4NURBS* G4Sphere::CreateNURBS () const
{
    return new G4NURBSbox (fRmax, fRmax, fRmax);       // Box for now!!!
}


// ******************************  End of G4Sphere.cc  ****************************************
