// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Torus.cc,v 1.4 1999-12-15 14:50:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Torus
//
// Implementation
//
// 30.10.96 V. Grichine First implementation with G4Tubs elements in Fs
// 09.10.98 V. Grichine modifications in Distance ToOut(p,v,...)
// 19.11.99 V. Grichine side = kNull in Distance ToOut(p,v,...)

#include "G4Torus.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBStube.hh"
#include "G4NURBScylinder.hh"
#include "G4NURBStubesector.hh"
#include "G4VisExtent.hh"


///////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4Torus::G4Torus(const G4String &pName,
	       G4double pRmin,
	       G4double pRmax,
	       G4double pRtor,
	       G4double pSPhi,
	       G4double pDPhi)
    : G4CSGSolid(pName)
{
    SetAllParameters(pRmin, pRmax, pRtor, pSPhi, pDPhi);
}

void
G4Torus::SetAllParameters(
	       G4double pRmin,
	       G4double pRmax,
	       G4double pRtor,
	       G4double pSPhi,
	       G4double pDPhi)
{
// Check swept radius
    if (pRtor>=pRmax)
	{
	   fRtor=pRtor;
	}
    else
    {
      G4Exception("Error in G4Torus::SetAllParameters - invalid swept radius");
    }

// Check radii

    if (pRmin<pRmax&&pRmin>=0)
    {
	   fRmin=pRmin; fRmax=pRmax;
    }
    else
    {
      G4Exception("Error in G4Torus::SetAllParameters - invalid radii");
    }

// Check angles
    if (pDPhi>=2.0*M_PI)
	{
	    fDPhi=2*M_PI;
	}
    else
    {
       if (pDPhi>0)
       {
		    fDPhi = pDPhi;
       }
       else
       {
	  G4Exception("Error in G4Torus::SetAllParameters - invalid dphi");
       }
    }
	
// Ensure psphi in 0-2PI or -2PI-0 range if shape crosses 0

    fSPhi = pSPhi;

    if (fSPhi<0)
	{
	    fSPhi=2.0*M_PI-fmod(fabs(fSPhi),2.0*M_PI);
	}
    else
	{
	    fSPhi=fmod(fSPhi,2.0*M_PI);
	}

    if (fSPhi+fDPhi>2.0*M_PI)
	{
	    fSPhi-=2.0*M_PI;
	}
}

//////////////////////////////////////////////////////////////////////
//
// Destructor

G4Torus::~G4Torus()
{;}

//////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Torus::ComputeDimensions(G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}

///////////////////////////////////////////////////////////////////////////
//
// Test function for study of intersections of a ray (starting from p along
// v) with the torus

G4int  G4Torus::TorusRoots(       G4double Ri,
			    const G4ThreeVector& p,
			    const G4ThreeVector& v) const
{
   // Define roots  Si (generally real >=0) for intersection with
   // torus (Ri = fRmax or fRmin) of ray p +S*v . General equation is :
   // c[4]*S^4 + c[3]*S^3 +c[2]*S^2 + c[1]*S + c[0] = 0 .
   
   G4double c[5],s[4] ;
   G4int num, i, j ;
   G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
   G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;
   G4double Rtor2 = fRtor*fRtor, Ri2 = Ri*Ri ;
   
   c[4] = 1.0 ;
   c[3] = 4*pDotV ;
   c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Ri2 + 2*Rtor2*v.z()*v.z()) ;
   c[1] = 4*(pDotV*(pRad2-Rtor2-Ri2) + 2*Rtor2*p.z()*v.z()) ;
   c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Ri2) 
          + 4*Rtor2*p.z()*p.z() + (Rtor2-Ri2)*(Rtor2-Ri2) ;
   
   num = SolveBiQuadratic(c,s) ;
   
   if(num)
   {
     for(i=0;i<num;i++)   // leave only >=0 roots
     {
	if(s[i]<0)
	{
	   for(j=i+1;j<num;j++) s[j-1] = s[j] ;
	   i-- ;
	   num-- ;
	}
     }
     if(num)
     {
        for(i=0;i<num;i++)
        {
           G4cout<<i<<" Root = "<<s[i]<<G4endl ; 
        }
     }
     else G4cout<<"All real roots are negative"<<G4endl ;
   }
   else G4cout<<"No real roots for intesection with torus"<<G4endl;
      
   return num ;      
}

/////////////////////////////////////////////////////////////////////////
//
// Auxiliary method for solving (in real numbers) biquadratic equation
// Algorithm based on : Graphics Gems I by Jochen Schwartz

G4int G4Torus::SolveBiQuadratic(double c[], double s[]  ) const
{
    G4double  coeffs[ 4 ];
    G4double  z, u, v, sub;
    G4double  A, B, C, D;
    G4double  A2, p, q, r;
    G4int     i,j, num;

    // normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 

    A = c[ 3 ];  // c[ 4 ]; since always c[4]==1 !
    B = c[ 2 ];  // c[ 4 ];
    C = c[ 1 ];  // c[ 4 ];
    D = c[ 0 ];  // c[ 4 ];

    //  substitute x = y - A/4 to eliminate cubic term:
    // y^4 + py^2 + qy + r = 0 

    A2 = A*A;
    p = - 0.375*A2 + B;   
    q = 0.125*A2*A - 0.5*A*B + C;
    r = - 3.0/256*A2*A2 + 1.0/16*A2*B - 0.25*A*C + D;

    // y^4 + py^2 + r = 0 and z=y^2 so y = +-sqrt(z1) and y = +-sqrt(z2)
   
    if(q==0) 
    {
        coeffs[ 0 ] = r;
	coeffs[ 1 ] = p;
	coeffs[ 2 ] = 1;
	num = SolveQuadratic(coeffs, s) ;

        if(num)
	{
          if(num==2)
	  {
	    if(s[0]>=0)
            {
	       if(s[0]==0) // Three roots and one of them == 0
	       {
	         s[2] = sqrt(s[1]) ;
	         s[1] = s[0] ;
		 s[0] = -s[2] ;
	         num++ ;
	       }
	       else        // Four roots
	       {
	         s[2] = sqrt(s[0]) ;
	         s[3] = sqrt(s[1]) ;
	         s[0] = -s[3] ;
	         s[1] = -s[2] ;
	         num +=2 ;
	       }
            }
	    else if(s[1]>=0)
	    {
	       if(s[1]==0)   // One root == 0
	       {
		  s[0] = 0 ;
		  num--;
	       }
	       else          // Two roots
	       {
	       s[0] = -sqrt(s[1]) ;
	       s[1] = -s[0] ;
	       }
	    }
	    else return num = 0 ; // Both Quadratic roots are negative
	  }
	  else    // num = 1 two equal roots from SolveQuadratic
	  {
	    if(s[0]>=0)
            {
               if(s[0]==0) ; 
	       else
	       {
	         s[1] = sqrt(s[0]) ;
	         s[0] = -s[1] ;
	         num +=1 ;
	       }
	    }
	    else return num = 0 ;
	  }
	}
	else return num ;
    }
    else if (r==0)	   // no absolute term: y(y^3 + py + q) = 0 
    {
	coeffs[ 0 ] = q;
	coeffs[ 1 ] = p;
	coeffs[ 2 ] = 0;
	coeffs[ 3 ] = 1;
	num = SolveCubic(coeffs, s);
	s[ num++ ] = 0;

	for(j=1;j<num;j++) // picksort of roots in ascending order
	{
	   sub = s[j] ;
	   i=j-1 ;
	   while(i>=0 && s[i]>sub)
	   {
	      s[i+1] = s[i--] ;
	   }
	   s[i+1] = sub ;
	}
    }
    else
    {
	// solve the resolvent cubic ... 

	coeffs[ 0 ] = 0.5*r*p - 0.125*q*q;
	coeffs[ 1 ] = - r;
	coeffs[ 2 ] = - 0.5*p;
	coeffs[ 3 ] = 1;

	num = SolveCubic(coeffs, s);

	// ... and take the one real solution ... 

	z = s[ 0 ];

	// ... to Build two quadratic equations 

	u = z * z - r;
	v = 2 * z - p;

	if (u==0)        u = 0 ;
	else if (u > 0)  u = sqrt(u) ;
	else             return 0 ;

	if (v==0)        v = 0 ;
	else if (v > 0)  v = sqrt(v);
	else             return 0 ;

	coeffs[ 0 ] = z - u;
	coeffs[ 1 ] = q < 0 ? -v : v;
	coeffs[ 2 ] = 1;

	num = SolveQuadratic(coeffs, s);

	coeffs[ 0 ]= z + u;
	coeffs[ 1 ] = q < 0 ? v : -v;
	coeffs[ 2 ] = 1;

	num += SolveQuadratic(coeffs, s + num);
    }

    // resubstitute 

    sub = 1.0/4 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;

    return num;
}

/////////////////////////////////////////////////////////////////////////////
//
// Auxiliary method for solving of cubic equation in real numbers
// From Graphics Gems I bu Jochen Schwartz

G4int G4Torus::SolveCubic(double c[], double s[]  ) const
{
    G4int     i, num;
    G4double  sub;
    G4double  A, B, C;
    G4double  A2, p, q;
    G4double  p3, D;

    // normal form: x^3 + Ax^2 + Bx + C = 0 

    A = c[ 2 ];           // c[ 3 ]; since always c[3]==1 !
    B = c[ 1 ];           // c[ 3 ];
    C = c[ 0 ];           // c[ 3 ];

    //  substitute x = y - A/3 to eliminate quadric term:
    //	x^3 +px + q = 0 

    A2 = A*A;
    p = 1.0/3*(- 1.0/3*A2 + B);
    q = 1.0/2*(2.0/27*A*A2 - 1.0/3*A*B + C);

    // use Cardano's formula 

    p3 = p*p*p;
    D = q*q + p3;

    if (D==0)
    {
	if (q==0) // one triple solution 
	{
	    s[ 0 ] = 0;
	    num = 1;
	}
	else // one single and one double solution 
	{
	    G4double u = cbrt(-q);
	    s[ 0 ] = 2 * u;
	    s[ 1 ] = - u;
	    num = 2;
	}
    }
    else if (D < 0) // Casus irreducibilis: three real solutions
    {
	G4double phi = 1.0/3 * acos(-q / sqrt(-p3));
	G4double t = 2 * sqrt(-p);

	s[ 0 ] =   t * cos(phi);
	s[ 1 ] = - t * cos(phi + M_PI / 3);
	s[ 2 ] = - t * cos(phi - M_PI / 3);
	num = 3;
    }
    else // one real solution 
    {
	G4double sqrt_D = sqrt(D);
	G4double u = cbrt(sqrt_D - q);
	G4double v = - cbrt(sqrt_D + q);

	s[ 0 ] = u + v;
	num = 1;
    }

    // resubstitute 

    sub = 1.0/3 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;

    return num;
}

///////////////////////////////////////////////////////////////////////////
//
// Auxiliary method for solving quadratic equations in real numbers
// From Graphics Gems I by Jochen Schwartz

G4int G4Torus::SolveQuadratic(double c[], double s[] ) const
{
    G4double p, q, D;

    // normal form: x^2 + px + q = 0 

    p = c[ 1 ]/2 ;             // * c[ 2 ]); since always c[2]==1
    q = c[ 0 ] ;               // c[ 2 ];

    D = p * p - q;

    if (D==0)
    {
	s[ 0 ] = - p;  // Generally we have two equal roots ?!
	return 1;      // But consider them as one for geometry
    }
    else if (D > 0)
    {
	G4double sqrt_D = sqrt(D);

	s[ 0 ] = - p - sqrt_D ;  // in ascending order !
	s[ 1 ] = - p + sqrt_D ;
	return 2;
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Torus::CalculateExtent(const EAxis pAxis,
			      const G4VoxelLimits& pVoxelLimit,
			      const G4AffineTransform& pTransform,
			      G4double& pMin, G4double& pMax) const
{
    if (!pTransform.IsRotated()&&fDPhi==2.0*M_PI&&fRmin==0)
	{
// Special case handling for unrotated solid torus
// Compute x/y/z mins and maxs for bounding box respecting limits,
// with early returns if outside limits. Then switch() on pAxis,
// and compute exact x and y limit for x/y case
	    
	    G4double xoffset,xMin,xMax;
	    G4double yoffset,yMin,yMax;
	    G4double zoffset,zMin,zMax;

	    G4double diff1,diff2,maxDiff,newMin,newMax;
	    G4double xoff1,xoff2,yoff1,yoff2;

	    xoffset=pTransform.NetTranslation().x();
	    xMin=xoffset-fRmax-fRtor;
	    xMax=xoffset+fRmax+fRtor;
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
	    yMin=yoffset-fRmax-fRtor;
	    yMax=yoffset+fRmax+fRtor;
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
	    zMin=zoffset-fRmax;
	    zMax=zoffset+fRmax;
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

// Known to cut cylinder
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
    else
	{
	    G4int i,noEntries,noBetweenSections4;
	    G4bool existsAfterClip=false;

// Calculate rotated vertex coordinates
	    G4ThreeVectorList *vertices;
            G4int noPolygonVertices ;  // will be 4 
	    vertices=CreateRotatedVertices(pTransform,noPolygonVertices);

	    pMin=+kInfinity;
	    pMax=-kInfinity;

	    noEntries=vertices->entries();
	    noBetweenSections4=noEntries-noPolygonVertices;
	    
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

////////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Torus::Inside(const G4ThreeVector& p) const
{
    G4double r2,pt2,pPhi,tolRMin,tolRMax;
    EInside in=kOutside;
                                              // General precals
            r2=p.x()*p.x()+p.y()*p.y();
	    pt2 = r2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*sqrt(r2) ;
	    if (fRmin) tolRMin=fRmin+kRadTolerance*0.5;
	    else tolRMin=0;
	    tolRMax=fRmax-kRadTolerance*0.5;
	    
	    if (pt2>=tolRMin*tolRMin && pt2<=tolRMax*tolRMax)
		{
		    if (fDPhi==2*M_PI||pt2==0)  // on torus swept axis
			{
			    in=kInside;
			}
		    else
			{
// Try inner tolerant phi boundaries (=>inside)
// if not inside, try outer tolerant phi boundaries
			    pPhi=atan2(p.y(),p.x());
			    if (pPhi<0) pPhi+=2*M_PI; // 0<=pPhi<2*M_PI
			    if (fSPhi>=0)
				{
				    if (pPhi>=fSPhi+kAngTolerance*0.5 &&
					pPhi<=fSPhi+fDPhi-kAngTolerance*0.5)
					{
					    in=kInside;
					}
				    else if (pPhi>=fSPhi-kAngTolerance*0.5 &&
					     pPhi<=fSPhi+fDPhi+kAngTolerance*0.5)
					{
					    in=kSurface;
					}
				}
			    else
				{
				    if (pPhi<fSPhi+2*M_PI) pPhi+=2*M_PI;
				    if (pPhi>=fSPhi+2*M_PI+kAngTolerance*0.5 &&
					pPhi<=fSPhi+fDPhi+2*M_PI-kAngTolerance*0.5)
					{
					    in=kInside;
					}
				    else if (pPhi>=fSPhi+2*M_PI-kAngTolerance*0.5 &&
					     pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5)
					{
					    in=kSurface;
					}
				}
			    
			    
			}
		}
	    else
		{
// Try generous boundaries
		    tolRMin=fRmin-kRadTolerance*0.5;
		    tolRMax=fRmax+kRadTolerance*0.5;
		    if (tolRMin<0) tolRMin=0;
		    if (pt2>=tolRMin*tolRMin && pt2 <= tolRMax*tolRMax)
			{
			    if (fDPhi==2*M_PI||pt2==0)
				{
// Continuous in phi or on z-axis
				    in=kSurface;
				}
			    else
				{
// Try outer tolerant phi boundaries only
				    pPhi=atan2(p.y(),p.x());
				    if (pPhi<0) pPhi+=2*M_PI; // 0<=pPhi<2*M_PI
				    if (fSPhi>=0)
					{
					    if (pPhi>=fSPhi-kAngTolerance*0.5 &&
						pPhi<=fSPhi+fDPhi+kAngTolerance*0.5)
						{
						    in=kSurface;
						}
					}
				    else
					{
					    if (pPhi<fSPhi+2*M_PI) pPhi+=2*M_PI;
					    if (pPhi>=fSPhi+2*M_PI-kAngTolerance*0.5 &&
						pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5)
						{
						    in=kSurface;
						}
					}
				    

				}
			}
		}
    return in;
}

/////////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Torus::SurfaceNormal( const G4ThreeVector& p) const
{
    ENorm side;
    G4ThreeVector norm;
    G4double rho2,rho,pt2,pt,phi;
    G4double distRMin,distRMax,distSPhi,distEPhi,distMin;

    rho2 = p.x()*p.x() + p.y()*p.y();
    rho = sqrt(rho2) ;
    pt2 = fabs(rho2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*rho) ;
    pt = sqrt(pt2) ;

    distRMax=fabs(pt-fRmax);

// First minimum
    if(fRmin)
    {
       distRMin=fabs(pt-fRmin);
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
	    norm=G4ThreeVector(-p.x()*(1-fRtor/rho)/pt,
			       -p.y()*(1-fRtor/rho)/pt,
			       -p.z()/pt);
	    break;
	case kNRMax:			// Outer radius
	    norm=G4ThreeVector(p.x()*(1-fRtor/rho)/pt,
			       p.y()*(1-fRtor/rho)/pt,
			       p.z()/pt);
	    break;
	case kNSPhi:
	    norm=G4ThreeVector(sin(fSPhi),-cos(fSPhi),0);
	    break;
	case kNEPhi:
	    norm=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
	    break;
	default:
	    G4Exception("Logic error in G4Torus::SurfaceNormal");
	    break;
	    
	} // end case

    return norm;
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

G4double G4Torus::DistanceToIn(const G4ThreeVector& p,
			       const G4ThreeVector& v) const
{
    G4double snxt=kInfinity, sphi=kInfinity;// snxt = default return value
    G4double c[5], s[4] ;

// Precalculated trig for phi intersections - used by r,z intersections to
//                                            check validity

    G4bool seg;				// true if segmented
    G4double hDPhi,hDPhiOT,hDPhiIT,cosHDPhiOT,cosHDPhiIT;
					// half dphi + outer tolerance
    G4double cPhi,sinCPhi,cosCPhi;	// central phi

    G4double tolORMin2,tolIRMin2;	// `generous' radii squared
    G4double tolORMax2,tolIRMax2 ;

    G4double Dist,xi,yi,zi,rhoi2,it2,inum,cosPsi; // Intersection point variables


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

// Calculate tolerant rmin and rmax
    if (fRmin>kRadTolerance)
	{
	    tolORMin2=(fRmin-0.5*kRadTolerance)*(fRmin-0.5*kRadTolerance);
	    tolIRMin2=(fRmin+0.5*kRadTolerance)*(fRmin+0.5*kRadTolerance);
	}
    else
	{
	    tolORMin2=0;
	    tolIRMin2=0;
	}
    tolORMax2=(fRmax+0.5*kRadTolerance)*(fRmax+0.5*kRadTolerance);
    tolIRMax2=(fRmax-kRadTolerance*0.5)*(fRmax-kRadTolerance*0.5);


//
// Intersection with Rmax (possible return) and Rmin (must also check phi)
//
   G4int    i,j,num ;
   G4double Rtor2=fRtor*fRtor, Rmax2=fRmax*fRmax, Rmin2=fRmin*fRmin ;
   G4double rho2 = p.x()*p.x()+p.y()*p.y();
   G4double rho = sqrt(rho2) ;
   G4double pt2 = fabs(rho2+p.z()*p.z() +Rtor2 - 2*fRtor*rho) ;
   //   G4double pt = sqrt(pt2) ;
   G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
   G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;
   G4double vDotNmax = pDotV -fRtor*(v.x()*p.x() + v.y()*p.y())/rho ;

// Inside outer radius :
// check not inside, and heading through tubs (-> 0 to in)
     if(pt2<=tolORMax2 && pt2>=tolIRMin2 && vDotNmax<0)
     {
	 if (seg)
   	 {
            inum = p.x()*cosCPhi+p.y()*sinCPhi;
	    cosPsi = inum/rho;
	    if (cosPsi>=cosHDPhiIT)
	    {
	       return snxt = 0;
	    }
	 }
	 else
	 {
	     return snxt = 0;
	 }
     }
     else         // intersection with Rmax torus
     {    
       c[4] = 1.0 ;
       c[3] = 4*pDotV ;
       c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmax2 + 2*Rtor2*v.z()*v.z()) ;
       c[1] = 4*(pDotV*(pRad2-Rtor2-Rmax2) + 2*Rtor2*p.z()*v.z()) ;
       c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmax2) 
              + 4*Rtor2*p.z()*p.z() + (Rtor2-Rmax2)*(Rtor2-Rmax2) ;
   
       num = SolveBiQuadratic(c,s) ;
   
       if(num)
       {
          for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots   P?!
          {
     	    if(s[i]<kRadTolerance*0.5)
            {
	      for(j=i+1;j<num;j++) s[j-1] = s[j] ;
	      i-- ;
	      num-- ;
	    }
          }
	  if(num)
	  {
	     for(i=0;i<num;i++)
	     {
	        if (seg)  // intersection point must have proper Phi
   	        {
		   xi=p.x()+s[i]*v.x();
		   yi=p.y()+s[i]*v.y();
 		   rhoi2=xi*xi+yi*yi;
                   inum = xi*cosCPhi + yi*sinCPhi;
	           cosPsi = inum/sqrt(rhoi2);
	           if (cosPsi>=cosHDPhiIT)
	           {
	               snxt = s[i] ;
		       break ;
	           }
	        }
	        else
	        {
	               snxt = s[i] ;
		       break ;
	        }
	     }
	  }
       }
     }
        
  if (fRmin)	// Possible Rmin intersection
  {
// Inside relative to inner radius :
// check not inside, and heading through tubs (-> 0 to in)
     if(pt2>=tolORMin2 && pt2<=tolIRMax2 && vDotNmax>0)
     {
	 if (seg)
   	 {
            inum = p.x()*cosCPhi+p.y()*sinCPhi;
	    cosPsi = inum/rho;
	    if (cosPsi>=cosHDPhiIT)
	    {
	       return snxt = 0;
	    }
	 }
	 else
	 {
	       return snxt = 0;
	 }
     }
     else              // intersection with Rmin torus
     {               
             c[4] = 1.0 ;
             c[3] = 4*pDotV ;
             c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmin2
	            + 2*Rtor2*v.z()*v.z()) ;
             c[1] = 4*(pDotV*(pRad2-Rtor2-Rmin2) + 2*Rtor2*p.z()*v.z()) ;
             c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmin2) 
                    + 4*Rtor2*p.z()*p.z() + (Rtor2-Rmin2)*(Rtor2-Rmin2) ;
   
             num = SolveBiQuadratic(c,s) ;
   
             if(num)
             {
               for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots   P?!
               {
     	         if(s[i]<kRadTolerance*0.5)
                 {
	            for(j=i+1;j<num;j++) s[j-1] = s[j] ;
	            i-- ;
	            num-- ;
	         }
               }
	       if(num)
	       {
	          for(i=0;i<num;i++)
	          {
	              if (seg)    // intersection point must have proper Phi
   	              {
		          xi=p.x()+s[i]*v.x();
		          yi=p.y()+s[i]*v.y();
 		          rhoi2=xi*xi+yi*yi;
                          inum = xi*cosCPhi + yi*sinCPhi;
	                  cosPsi = inum/sqrt(rhoi2);
	                  if (cosPsi>=cosHDPhiIT && s[i]<snxt)
	                  {
	                     snxt = s[i] ;
		             break ;
	                  }
	              }
	              else if(s[i]<snxt)
	              {
	                      snxt = s[i] ;
		              break ;
	              }
	           }
	       }
             }
        }
    }      // if(Rmin)
		
//
// Phi segment intersection
//
// o Tolerant of points inside phi planes by up to kCarTolerance*0.5
//
// o NOTE: Large duplication of code between sphi & ephi checks
//         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
//            intersection check <=0 -> >=0
//         -> use some form of loop Construct ?
//
    if (seg)
    {
                                      // First phi surface (`S'tarting phi)
    sinSPhi=sin(fSPhi);
    cosSPhi=cos(fSPhi);
    Comp=v.x()*sinSPhi-v.y()*cosSPhi;  // Compnent in outwards normal dirn
	                	
	    if (Comp<0)
		{
		    Dist=(p.y()*cosSPhi-p.x()*sinSPhi);
		    if (Dist<kCarTolerance*0.5)
		    {
			 sphi=Dist/Comp;
			 if (sphi<snxt)
			 {
			    if (sphi<0)
			    {
			       sphi=0;
			    }
                	    xi=p.x()+sphi*v.x();
			    yi=p.y()+sphi*v.y();
			    zi=p.z()+sphi*v.z();
 			    rhoi2=xi*xi+yi*yi;
                            it2 = fabs(rhoi2+zi*zi +Rtor2 - 2*fRtor*sqrt(rhoi2)) ;
 			    if (it2>=tolORMin2 && it2<=tolORMax2)
 			    {
//  r intersection is good - check intersecting with correct half-plane

			       if ((yi*cosCPhi-xi*sinCPhi)<=0)	snxt=sphi;
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
		    if (Dist<kCarTolerance*0.5)
		    {
		      sphi=Dist/Comp;
		      if (sphi<snxt)
		      {
			 if (sphi<0)
			 {
			    sphi=0;
			 }
                	 xi=p.x()+sphi*v.x();
			 yi=p.y()+sphi*v.y();
			 zi=p.z()+sphi*v.z();
 			 rhoi2=xi*xi+yi*yi;
                         it2 = fabs(rhoi2+zi*zi +Rtor2 - 2*fRtor*sqrt(rhoi2)) ;

 			 if (it2>=tolORMin2 && it2<=tolORMax2)
 			 {
// z and r intersections good - check intersecting with correct half-plane

			    if ((yi*cosCPhi-xi*sinCPhi)>=0)	snxt=sphi;
 			 }    
		       }
		     }
		}
	}         // if(seg)   
	

    return snxt;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Torus::DistanceToIn(const G4ThreeVector& p) const
{
    G4double safe,safe1,safe2;
    G4double phiC,cosPhiC,sinPhiC,safePhi,ePhi,cosPsi;
    G4double rho2,rho,pt2,pt ;
    
    rho2 = p.x()*p.x()+p.y()*p.y();
    rho  = sqrt(rho2) ;
    pt2  = fabs(rho2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*rho) ;
    pt   = sqrt(pt2) ;

    safe1=fRmin-pt;
    safe2=pt-fRmax;

    if (safe1>safe2) safe=safe1;
    else safe=safe2;

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

///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4Torus::DistanceToOut(const G4ThreeVector& p,
				const G4ThreeVector& v,
			        const G4bool calcNorm,
			        G4bool *validNorm,
				G4ThreeVector  *n    ) const
{
    ESide side = kNull, sidephi ;
    G4double snxt=kInfinity, sphi,c[5],s[4];

// Vars for phi intersection

    G4double sinSPhi,cosSPhi,ePhi,sinEPhi,cosEPhi;
    G4double cPhi,sinCPhi,cosCPhi;
    G4double pDistS,compS,pDistE,compE,sphi2,xi,yi,zi,vphi;

// Radial Intersections Defenitions & General Precals
    
   // Define roots  Si (generally real >=0) for intersection with
   // torus (Ri = fRmax or fRmin) of ray p +S*v . General equation is :
   // c[4]*S^4 + c[3]*S^3 +c[2]*S^2 + c[1]*S + c[0] = 0 .
   
   G4int    i,j,num ;
   G4double Rtor2=fRtor*fRtor, Rmax2=fRmax*fRmax, Rmin2=fRmin*fRmin ;
   G4double rho2 = p.x()*p.x()+p.y()*p.y();
   G4double rho = sqrt(rho2) ;
   G4double pt2 = fabs(rho2+p.z()*p.z() + Rtor2 - 2*fRtor*rho) ;
   G4double pt = sqrt(pt2) ;
   G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
   G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;
   
   G4double tolRMax=fRmax-kRadTolerance*0.5;
   
   G4double vDotNmax = pDotV -fRtor*(v.x()*p.x() + v.y()*p.y())/rho ;
   G4double pDotxyNmax = (1-fRtor/rho) ;

     
     if(pt2>tolRMax*tolRMax && vDotNmax>=0)
     {
// On tolerant boundary & heading outwards (or perpendicular to) outer
// radial surface -> leaving immediately with *n for really convex part only
	 if (calcNorm && pDotxyNmax>=-kRadTolerance) 
	 {
	    *n=G4ThreeVector(p.x()*(1-fRtor/rho)/pt,
	                     p.y()*(1-fRtor/rho)/pt,
	                     p.z()/pt);
	    *validNorm=true;
	 }
         return snxt=0; // Leaving by Rmax immediately
     }
     else
     {     // intersection with Rmax torus
       c[4] = 1.0 ;
       c[3] = 4*pDotV ;
       c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmax2 + 2*Rtor2*v.z()*v.z()) ;
       c[1] = 4*(pDotV*(pRad2-Rtor2-Rmax2) + 2*Rtor2*p.z()*v.z()) ;
       c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmax2) 
              + 4*Rtor2*p.z()*p.z() + (Rtor2-Rmax2)*(Rtor2-Rmax2) ;
   
       num = SolveBiQuadratic(c,s) ;
   
       if(num)
       {
          for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots
          {
     	    if(s[i]<kRadTolerance*0.5)
            {
	      for(j=i+1;j<num;j++) s[j-1] = s[j] ;
	      i-- ;
	      num-- ;
	    }
          }
	  if(num)
	  {
	     snxt = s[0] ;
	     side = kRMax ;
	  }
       }
        	// Possible Rmin intersection
       if (fRmin)
       {
           G4double tolRMin=fRmin+kRadTolerance*0.5;
// Leaving via Rmin
// NOTE: SHould use rho-rmin>kRadTolerance*0.5 - avoid sqrt for efficiency
	   if (pt2<tolRMin*tolRMin && vDotNmax<0)
    	   {
	      if (calcNorm)
	      {
	        *validNorm=false;      // Concave surface of the torus
	      }
              return snxt=0;           // Leaving by Rmin immediately
	   }
           else
	   {                // intersection with Rmin torus
             c[4] = 1.0 ;
             c[3] = 4*pDotV ;
             c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmin2
	            + 2*Rtor2*v.z()*v.z()) ;
             c[1] = 4*(pDotV*(pRad2-Rtor2-Rmin2) + 2*Rtor2*p.z()*v.z()) ;
             c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmin2) 
                    + 4*Rtor2*p.z()*p.z() + (Rtor2-Rmin2)*(Rtor2-Rmin2) ;
   
             num = SolveBiQuadratic(c,s) ;
   
             if(num)
             {
               for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots
               {
     	         if(s[i]<kRadTolerance*0.5)
                 {
	            for(j=i+1;j<num;j++) s[j-1] = s[j] ;
	            i-- ;
	            num-- ;
	         }
               }
	       if(num && s[0]<snxt)
	       {
	          snxt = s[0] ;
	          side = kRMin ;
	       }
             }
	   }
        }      // if(Rmin)
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
						sphi=kInfinity;
					    else
						{
						    sidephi=kSPhi;
						   if (pDistS>-kCarTolerance*0.5)
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
						  if (pDistE<=-kCarTolerance*0.5)
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
							    if (pDistE>-kCarTolerance*0.5)
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
							    if (pDistS>-kCarTolerance*0.5)
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
    G4double rhoi2,rhoi,it2,it,iDotxyNmax ;

    if (calcNorm)
	{
	    switch(side)
		{
		case kRMax:               // n is unit vector 
		    xi=p.x()+snxt*v.x();
		    yi=p.y()+snxt*v.y();
		    zi=p.z()+snxt*v.z();
                    rhoi2 = xi*xi+yi*yi;
                    rhoi = sqrt(rhoi2) ;
                    it2 = fabs(rhoi2+zi*zi +fRtor*fRtor - 2*fRtor*rhoi) ;
                    it = sqrt(it2) ;
                    iDotxyNmax = (1-fRtor/rhoi) ;
		    if(iDotxyNmax>=-kRadTolerance)
		    {                       // really convex part of Rmax
	               *n=G4ThreeVector(xi*(1-fRtor/rhoi)/it,
	                                yi*(1-fRtor/rhoi)/it,
	                                zi/it);
		       *validNorm=true;
		    }
		    else
		    {
		       *validNorm=false; // concave-convex part of Rmax
		    }
		    break;
		case kRMin:
		    *validNorm=false;	// Rmin is concave or concave-convex
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
		default:
		    G4Exception("Invalid enum in G4Torus::DistanceToOut");
		    break;
		}
	}

    return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance (<=actual) to closest surface of shape from inside

G4double G4Torus::DistanceToOut(const G4ThreeVector& p) const
{
    G4double safe,safeR1,safeR2;
    G4double rho2,rho,pt2,pt ;
    G4double safePhi,phiC,cosPhiC,sinPhiC,ePhi;
    rho2=p.x()*p.x()+p.y()*p.y();
    rho=sqrt(rho2);
    pt2 = fabs(rho2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*rho) ;
    pt = sqrt(pt2) ;

    if (fRmin)
	{
	    safeR1=pt-fRmin;
	    safeR2=fRmax-pt;
	    if (safeR1<safeR2)
		{
		    safe=safeR1;
		}
	    else
		{
		    safe=safeR2;
		}
	}
    else
	{
	    safe=fRmax-pt;
	}
    

// Check if phi divided, Calc distances closest phi plane
    if (fDPhi<2.0*M_PI)
	{
// Above/below central phi of Torus?
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
   G4Torus::CreateRotatedVertices(const G4AffineTransform& pTransform,
				  G4int& noPolygonVertices) const
{
    G4ThreeVectorList *vertices;
    G4ThreeVector vertex0,vertex1,vertex2,vertex3;
    G4double meshAngle,meshRMax,crossAngle,cosCrossAngle,sinCrossAngle,sAngle;
    G4double rMaxX,rMaxY,rMinX,rMinY;
    G4int crossSection,noCrossSections;

// Compute no of cross-sections necessary to mesh tube

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
    meshRMax=(fRtor+fRmax)/cos(meshAngle*0.5);

// If complete in phi, set start angle such that mesh will be at fRmax
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

		    rMaxX=meshRMax*cosCrossAngle;
		    rMaxY=meshRMax*sinCrossAngle;
		    rMinX=(fRtor-fRmax)*cosCrossAngle;
		    rMinY=(fRtor-fRmax)*sinCrossAngle;
		    vertex0=G4ThreeVector(rMinX,rMinY,-fRmax);
		    vertex1=G4ThreeVector(rMaxX,rMaxY,-fRmax);
		    vertex2=G4ThreeVector(rMaxX,rMaxY,+fRmax);
		    vertex3=G4ThreeVector(rMinX,rMinY,+fRmax);

		    vertices->insert(pTransform.TransformPoint(vertex0));
		    vertices->insert(pTransform.TransformPoint(vertex1));
		    vertices->insert(pTransform.TransformPoint(vertex2));
		    vertices->insert(pTransform.TransformPoint(vertex3));
		}
                noPolygonVertices = 4 ;
	}
    else
	{
	    G4Exception("G4Torus::CreateRotatedVertices Out of memory - Cannot alloc vertices");
	}
    return vertices;
}

///////////////////////////////////////////////////////////////////////
//
// No implementation for Visualisation Functions

void G4Torus::DescribeYourselfTo (G4VGraphicsScene& scene) const {
  scene.AddThis (*this);
}

G4VisExtent G4Torus::GetExtent() const {
  // Define the sides of the box into which the G4Torus instance would fit.
  return G4VisExtent (-fRtor-fRmax, fRtor+fRmax,
		      -fRtor-fRmax, fRtor+fRmax, -fRmax, fRmax);
}

G4Polyhedron* G4Torus::CreatePolyhedron () const {
  return new G4PolyhedronTorus (fRmin, fRmax, fRtor, fSPhi, fSPhi + fDPhi);
}

G4NURBS* G4Torus::CreateNURBS () const {
  G4NURBS* pNURBS;
  if (fRmin != 0) {
    if (fDPhi >= 2.0 * M_PI) {
      pNURBS = new G4NURBStube (fRmin, fRmax, fRtor);
    }
    else {
      pNURBS = new G4NURBStubesector (fRmin, fRmax, fRtor, fSPhi, fSPhi + fDPhi);
    }
  }
  else {
    if (fDPhi >= 2.0 * M_PI) {
      pNURBS = new G4NURBScylinder (fRmax, fRtor);
    }
    else {
      const G4double epsilon = 1.e-4; // Cylinder sector not yet available!
      pNURBS = new G4NURBStubesector (epsilon, fRmax, fRtor,
				      fSPhi, fSPhi + fDPhi);
    }
  }
  return pNURBS;
}

//
//
/////////////////////////////////////////////////////////////////////////////


