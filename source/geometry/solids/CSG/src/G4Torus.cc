// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Torus.cc,v 1.23 2001-01-29 13:12:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Torus
//
// Implementation
//
// 30.10.96 V.Grichine First implementation with G4Tubs elements in Fs
// 09.10.98 V.Grichine modifications in Distance ToOut(p,v,...)
// 19.11.99 V.Grichine side = kNull in Distance ToOut(p,v,...)
// 06.03.00 V.Grichine, modifications in Distance ToOut(p,v,...)
// 26.05.00 V.Grichine, new fuctions developed by O.Cremonesi were added
// 31.08.00 E.Medernach, numerical computation of roots with bounding volume technique
// 03.10.00 E.Medernach, SafeNewton added
// 11.01.01 E.Medernach, Use G4PolynomialSolver to find roots
//


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
#include "G4PolynomialSolver.hh"

// #define DEBUGTORUS 1

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
  if ( pRtor >= pRmax + kCarTolerance )      // Check swept radius
  {
    fRtor = pRtor ;
  }
  else
  {
    G4Exception("Error in G4Torus::SetAllParameters - invalid swept radius");
  }

// Check radii

  if (pRmin < pRmax - 2*kCarTolerance && pRmin >= 0 )
  {
    if (pRmin >= kCarTolerance) fRmin = pRmin ;
    else                        fRmin = 0.0   ;
 
    fRmax = pRmax ;
  }
  else
  {
    G4Exception("Error in G4Torus::SetAllParameters - invalid radii");
  }

// Check angles

  if ( pDPhi >= 2.0*M_PI )  fDPhi = 2*M_PI ;
  else
  {
    if (pDPhi > 0)   fDPhi = pDPhi ;
    else
    {
      G4Exception("Error in G4Torus::SetAllParameters - invalid dphi");
    }
  }
	
// Ensure psphi in 0-2PI or -2PI-0 range if shape crosses 0

  fSPhi = pSPhi;

  if (fSPhi < 0) fSPhi = 2.0*M_PI - fmod(fabs(fSPhi), 2.0*M_PI) ;

  else fSPhi = fmod(fSPhi, 2.0*M_PI) ;

  if (fSPhi+fDPhi > 2.0*M_PI) fSPhi -= 2.0*M_PI ;
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
 
   num = SolveBiQuadraticNew(c,s) ;
   
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
           G4cout<<i<<" new Root = "<<s[i]<<G4endl ; 
        }
     }
     else G4cout<<"All real new roots are negative"<<G4endl ;
   }
   else G4cout<<"No real new roots for intesection with torus"<<G4endl;

      
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
    else if (r == 0)	   // no absolute term: y(y^3 + py + q) = 0 
    {
	coeffs[ 0 ] = q ;
	coeffs[ 1 ] = p ;
	coeffs[ 2 ] = 0 ;
	coeffs[ 3 ] = 1 ;
	num = SolveCubic(coeffs, s) ;

	s[ num++ ] = 0;

	for(j=1;j<num;j++) // picksort of roots in ascending order
	{
	   sub = s[j] ;
	   i=j-1 ;
	   while( i >= 0 && s[i] > sub )
	   {
	      i-- ;    
	      s[i+1] = s[i] ;           // s[i--] ;
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

    if (D == 0)
    {
	if (q == 0) // one triple solution 
	{
	    s[ 0 ] = 0;
	    num = 1;
	}
	else // one single and one double solution 
	{
	    G4double u = pow(-q,1./3.);
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
	G4double u = pow(sqrt_D - q,1./3.);
	G4double v = - pow(sqrt_D + q,1./3.);

	s[ 0 ] = u + v;
	num = 1;
    }

    // resubstitute 

    sub = 1.0/3 * A;

    for (i = 0; i < num; ++i)
	s[ i ] -= sub;

    return num;
}

// ---------------------------------------------------------------------

G4int G4Torus::SolveBiQuadraticNew(double c[], double s[]  ) const
{
// From drte4 by McLareni; rewritten by O.Cremonesi

    G4double  coeffs[ 4 ];
    G4double  w1, w2, w3;
    G4double  sub;
    G4double  A, B, C, D;
    G4double  A2, p, q, r ;
    G4int     i,j, num;

    // normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 

    A = c[ 3 ];  // c[ 4 ]; since always c[4]==1 !
    B = c[ 2 ];  // c[ 4 ];
    C = c[ 1 ];  // c[ 4 ];
    D = c[ 0 ];  // c[ 4 ];

    if( B==0 && C==0 ) 
    {
      if( D==0 ) 
      {
 	s[0] = -A;
	s[1] = s[2] = s[3] = 0;
	return 4;
      }
    }
    else if( A==0 ) 
    {
      if( D>0 ) return 0;
      else 
      {
 	s[0] = sqrt( sqrt( -D ) );
 	s[1] = -s[0];
 	return 2;
      }
    }
    
    //  substitute x = y - A/4 to eliminate cubic term:
    // y^4 + py^2 + qy + r = 0 

    A2 = A*A;
    p = B - 3.0*A2/8.0;   
    q = C - 0.5*A*( B-A2/4.0 );
    r = D - (A*C-A2/4.0*(B-A2*3.0/16.0))/4.0;
    coeffs[ 0 ] = -q*q/64.;
    coeffs[ 1 ] = (p*p/4.0-r)/4.0;
    coeffs[ 2 ] = p/2.0;
    coeffs[ 3 ] = 1;
    
    G4double cubic_discr;
    num = SolveCubicNew(coeffs, s, cubic_discr);
    
    sub = A/4.0;
    num = 0;
    
    if( cubic_discr == 0 ) s[2] = s[1];
    
    if( cubic_discr <= 0 ) 
    {
      num = 4;
      G4double v[3];
      G4double vm1 = -1.0e99, vm2 ;
      for( i=0; i<3; i++ ) 
      {
        v[i] = fabs( s[i] ) ;
        if( v[i] > vm1 ) vm1 = v[i] ;
      }
      if( vm1 == v[0] ) 
      {
        i = 0;
        if( v[1] > v[2] ) vm2 = v[1];
        else vm2 = v[2];
      } 
      else if( vm1 == v[1] ) 
      {
        i = 1;
        if( v[0] > v[2] ) vm2 = v[0];
        else vm2 = v[2];
      }  
      else 
      {
        i = 2;
        if( v[0] > v[1] ) vm2 = v[0];
        else vm2 = v[1];
      }
      if( vm2 == v[0] )      j = 0 ;
      else if( vm2 == v[1] ) j = 1 ;
      else j = 2 ;

      w1 = sqrt( s[i] );
      w2 = sqrt( s[j] );
    } 
    else 
    {
     num = 2;
     w1 = w2 = sqrt( s[1] );
    }
    if( w1*w2 != 0. ) w3 = -q/( 8.0*w1*w2 ) ;
    else              w3 = 0.0 ;
    
    if( num == 4 ) 
    {
      s[0] =  w1 + w2 + w3 - sub ;
      s[1] = -w1 - w2 + w3 - sub ;
      s[2] = -w1 + w2 - w3 - sub ;
      s[3] =  w1 - w2 - w3 - sub ;
    }
    else if( num == 2 ) 
    {
      s[0] =  w1 + w2 + w3 - sub ;
      s[1] = -w1 - w2 + w3 - sub ;
    }     
    return num ;
}

// -------------------------------------------------------------------------

G4int G4Torus::SolveCubicNew(double c[], double s[], double& cubic_discr ) const
{
// From drte3 by McLareni; rewritten by O.Cremonesi
    const G4double eps = 1.e-6;
    const G4double delta = 1.e-15;
    G4int     i, j;
    G4double  sub;
    G4double  y[3];
    G4double  A, B, C;
    G4double  A2, p, q;
    G4double  h1,h2,h3;
    G4double  u,v;

    // normal form: x^3 + Ax^2 + Bx + C = 0 

    A = c[ 2 ];           // c[ 3 ]; since always c[3]==1 !
    B = c[ 1 ];           // c[ 3 ];
    C = c[ 0 ];           // c[ 3 ];

    if( B==0 && C==0 ) 
    {
      s[0] = -A;
      s[1] = s[2] = 0.;
      cubic_discr = 0.;
      return 3;
    }
    A2 = A*A;
    p = B - A2/3.0;
    q = ( A2*2.0/27.-B/3.0 )*A + C;
    cubic_discr = q*q/4.0 + p*p*p/27.0;
    sub = A/3.0;
    h1 = q/2.0;

    if( cubic_discr > delta ) 
    {
      h2 = sqrt( cubic_discr );
      u = -h1+h2;
      v = -h1-h2;
      if( u < 0 ) u = -pow(-u,1./3.);
      else u = pow(u,1./3.);
      if( v < 0 ) v = -pow(-v,1./3.);
      else v = pow(v,1./3.);
      s[0] = u+v-sub;
      s[1] = -(u+v)/2.0-sub;
      s[2] = fabs(u-v)*sqrt(3.0)/2.0;
      if( fabs(u) <= eps || fabs(v) <= eps ) 
      {
        y[0] = s[0] ;

	for( i=0; i<2; i++ ) 
	{
          y[i+1] = y[i] - (((y[i]+A)*y[i]+B)*y[i]+C)/((3.*y[i]+2.*A)*y[i]+B);
	}
	s[0] = y[2];
	return 1;
      }
    }
    else if( fabs(cubic_discr) <= delta ) 
    {
      cubic_discr = 0.;

      if( h1 < 0 ) u = pow(-h1,1./3.);
      else         u = -pow(h1,1./3.);

      s[0] =  u + u - sub ;
      s[1] = -u - sub ;
      s[2] = s[1] ;

      if( fabs(h1) <= eps ) 
      {
        y[0] = s[0];
 	for( i=0; i<2; i++ ) 
        {
	  h1 = (3.0*y[i]+2.*A)*y[i]+B;

	  if( fabs(h1) > delta ) y[i+1] = y[i]-(((y[i]+A)*y[i]+B)*y[i]+C)/h1;
	  else 
          {
	    s[0] = s[1] = s[2] = -A/3.;
	    return 3;
	  }
	}
	s[0] = y[2];
	s[1] = s[2] = -(A+s[0])/2.;
	return 3;
      } 
    }
    else 
    {
      h3 =fabs(p/3.);
      h3 = sqrt(h3*h3*h3);
      h2 = acos(-h1/h3)/3.;
      h1 = pow(h3,1./3.);
      u = h1*cos(h2);
      v = sqrt(3.)*h1*sin(h2);
      s[0] = u+u-sub;
      s[1] = -u-v-sub;
      s[2] = -u+v-sub;

      if( h3 <= eps || s[0] <=eps || s[1] <= eps || s[2] <= eps ) 
      {
        for( i=0; i<3; i++ ) 
        {
	  y[0] = s[i] ;

	  for( j=0; j<2; j++ )
	  {
	      y[j+1] = y[j]-(((y[j]+A)*y[j]+B)*y[j]+C)/((3.*y[j]+2.*A)*y[j]+B);
	  }
	  s[i] = y[2] ;
	}
      }
    }
    return 3;
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
  if (!pTransform.IsRotated() && fDPhi==2.0*M_PI && fRmin==0)
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

    xoffset = pTransform.NetTranslation().x();
    xMin    = xoffset - fRmax - fRtor ;
    xMax    = xoffset + fRmax + fRtor ;

    if (pVoxelLimit.IsXLimited())
    {
      if (xMin > pVoxelLimit.GetMaxXExtent()+kCarTolerance || 
          xMax < pVoxelLimit.GetMinXExtent()-kCarTolerance)     return false ;
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
      if (yMin > pVoxelLimit.GetMaxYExtent()+kCarTolerance || 
          yMax < pVoxelLimit.GetMinYExtent()-kCarTolerance) return false ;
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
      if (zMin > pVoxelLimit.GetMaxZExtent()+kCarTolerance || 
          zMax < pVoxelLimit.GetMinZExtent()-kCarTolerance   ) return false ;
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

	  pMin = xMin ;
	  pMax = xMax ;
	}
	else
        {
// Y limits don't cross max/min x => compute max delta x, hence new mins/maxs

	  diff1   = sqrt(fRmax*fRmax - yoff1*yoff1) ;
	  diff2   = sqrt(fRmax*fRmax - yoff2*yoff2) ;
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
			    
          pMin = yMin ;
	  pMax = yMax ;
	} 
	else
	{
// X limits don't cross max/min y => compute max delta y, hence new mins/maxs

          diff1   = sqrt(fRmax*fRmax - xoff1*xoff1) ;
	  diff2   = sqrt(fRmax*fRmax - xoff2*xoff2) ;
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

    noEntries          = vertices->entries() ;
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

////////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Torus::Inside(const G4ThreeVector& p) const
{
  G4double r2, pt2, pPhi, tolRMin, tolRMax ;

  EInside in = kOutside ;
                                              // General precals
  r2  = p.x()*p.x() + p.y()*p.y() ;
  pt2 = r2 + p.z()*p.z() + fRtor*fRtor - 2*fRtor*sqrt(r2) ;

  if (fRmin) tolRMin = fRmin + kRadTolerance*0.5 ;
  else       tolRMin = 0 ;

  tolRMax = fRmax - kRadTolerance*0.5;
	    
  if (pt2 >= tolRMin*tolRMin && pt2 <= tolRMax*tolRMax )
  {
    if ( fDPhi == 2*M_PI || pt2 == 0 )  // on torus swept axis
    {
      in = kInside ;
    }
    else
    {
// Try inner tolerant phi boundaries (=>inside)
// if not inside, try outer tolerant phi boundaries

      pPhi = atan2(p.y(),p.x()) ;

      if ( pPhi  <  0 ) pPhi += 2*M_PI ; // 0<=pPhi<2*M_PI
      if ( fSPhi >= 0 )
      {
        if ( pPhi >= fSPhi+kAngTolerance*0.5 &&
	     pPhi <= fSPhi+fDPhi-kAngTolerance*0.5 ) in = kInside ;

	else if ( pPhi >= fSPhi-kAngTolerance*0.5 &&
		  pPhi <= fSPhi+fDPhi+kAngTolerance*0.5 ) in = kSurface ;
      }
      else
      {
        if (pPhi < fSPhi+2*M_PI) pPhi += 2*M_PI ;

	if ( pPhi >= fSPhi+2*M_PI+kAngTolerance*0.5 &&
	     pPhi <= fSPhi+fDPhi+2*M_PI-kAngTolerance*0.5 ) in = kInside ;

	else if ( pPhi >= fSPhi+2*M_PI-kAngTolerance*0.5 &&
		  pPhi <= fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5) in = kSurface ;
      }			    
    }
  }
  else   // Try generous boundaries
  {
    tolRMin = fRmin - kRadTolerance*0.5 ;
    tolRMax = fRmax + kRadTolerance*0.5 ;

    if (tolRMin < 0 ) tolRMin = 0 ;

    if (pt2 >= tolRMin*tolRMin && pt2 <= tolRMax*tolRMax)
    {
      if (fDPhi == 2*M_PI || pt2 == 0 ) // Continuous in phi or on z-axis
      {
        in = kSurface ;
      }
      else // Try outer tolerant phi boundaries only
      {
        pPhi = atan2(p.y(),p.x()) ;

	if ( pPhi < 0 ) pPhi += 2*M_PI ; // 0<=pPhi<2*M_PI

	if (fSPhi >= 0 )
        {
	  if( pPhi >= fSPhi-kAngTolerance*0.5 &&
	      pPhi <= fSPhi+fDPhi+kAngTolerance*0.5) in = kSurface ;
	}
	else
        {
	  if (pPhi < fSPhi + 2*M_PI) pPhi += 2*M_PI ;

	  if ( pPhi >= fSPhi+2*M_PI-kAngTolerance*0.5 &&
	       pPhi <= fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5 ) in = kSurface ; 
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

G4ThreeVector G4Torus::SurfaceNormal( const G4ThreeVector& p) const
{
  ENorm side ;
  G4ThreeVector norm;
  G4double rho2,rho,pt2,pt,phi;
  G4double distRMin,distRMax,distSPhi,distEPhi,distMin;

  rho2 = p.x()*p.x() + p.y()*p.y();
  rho = sqrt(rho2) ;
  pt2 = fabs(rho2+p.z()*p.z() +fRtor*fRtor - 2*fRtor*rho) ;
  pt = sqrt(pt2) ;

  distRMax = fabs(pt - fRmax) ;


  if(fRmin)  // First minimum radius
  {
    distRMin = fabs(pt - fRmin) ;

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
  if (fDPhi < 2.0*M_PI && rho )
  {
    phi = atan2(p.y(),p.x()) ; // Protected against (0,0,z) (above rho !=0)

    if (phi < 0) phi += 2*M_PI ;

    if (fSPhi < 0 ) distSPhi = fabs(phi-(fSPhi+2.0*M_PI))*rho ;
    else            distSPhi = fabs(phi-fSPhi)*rho ;

    distEPhi = fabs(phi - fSPhi - fDPhi)*rho ;

    if (distSPhi < distEPhi) // Find new minimum
    {
      if (distSPhi<distMin) side = kNSPhi ;
    }
    else
    {
      if (distEPhi < distMin) side = kNEPhi ;
    }
  }	
  switch (side)
  {
    case kNRMin:			// Inner radius
      norm = G4ThreeVector( -p.x()*(1-fRtor/rho)/pt,
			    -p.y()*(1-fRtor/rho)/pt,
			    -p.z()/pt                 ) ;
      break ;
    case kNRMax:			// Outer radius
      norm = G4ThreeVector( p.x()*(1-fRtor/rho)/pt,
			    p.y()*(1-fRtor/rho)/pt,
			    p.z()/pt                  ) ;
      break;
    case kNSPhi:
      norm = G4ThreeVector(sin(fSPhi),-cos(fSPhi),0) ;
      break;
    case kNEPhi:
      norm = G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0) ;
      break;
    default:
      G4Exception("Logic error in G4Torus::SurfaceNormal");
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

G4double G4Torus::DistanceToIn(const G4ThreeVector& p,
			       const G4ThreeVector& v) const
{

  /*
    On voudrait arriver a cela:
    return SolveNumeric(p, v, true);
    Mais des problemes avec la tolerance sur la section Phi
    ne le permet pas pour le moment
  */

  /*
    Le tore mathematique peut etre vu comme une equation implicite
   */
  G4double snxt=kInfinity, sphi=kInfinity;// snxt = default return value

  G4double c[5], s[4] ;

// Precalculated trig for phi intersections - used by r,z intersections to
//                                            check validity

  G4bool seg;				// true if segmented
  G4double hDPhi,hDPhiOT,hDPhiIT,cosHDPhiOT=0.,cosHDPhiIT=0.;
					// half dphi + outer tolerance
  G4double cPhi,sinCPhi=0.,cosCPhi=0.;	// central phi

  G4double tolORMin2,tolIRMin2;	// `generous' radii squared
  G4double tolORMax2,tolIRMax2 ;

  G4double Dist,xi,yi,zi,rhoi2,it2,inum,cosPsi; // Intersection point variables


  G4double Comp;
  G4double cosSPhi,sinSPhi;		// Trig for phi start intersect
  G4double ePhi,cosEPhi,sinEPhi;	// for phi end intersect

#if DEBUGTORUS
	G4cout << "G4Torus::DistanceToIn    " << p << ", " << v << G4endl;
#endif

// Set phi divided flag and precalcs

  if ( fDPhi < 2.0*M_PI )
  {
    seg        = true ;
    hDPhi      = 0.5*fDPhi ;		// half delta phi
    cPhi       = fSPhi + hDPhi ;
    hDPhiOT    = hDPhi+0.5*kAngTolerance ;	// outers tol' half delta phi 
    hDPhiIT    = hDPhi - 0.5*kAngTolerance ;
    sinCPhi    = sin(cPhi) ;
    cosCPhi    = cos(cPhi) ;
    cosHDPhiOT = cos(hDPhiOT) ;
    cosHDPhiIT = cos(hDPhiIT) ;
  }
  else seg = false ;

  if (fRmin > kRadTolerance) // Calculate tolerant rmin and rmax
  {
    tolORMin2 = (fRmin - 0.5*kRadTolerance)*(fRmin - 0.5*kRadTolerance) ;
    tolIRMin2 = (fRmin + 0.5*kRadTolerance)*(fRmin + 0.5*kRadTolerance) ;
  }
  else
  {
    tolORMin2 = 0 ;
    tolIRMin2 = 0 ;
  }
  tolORMax2 = (fRmax + 0.5*kRadTolerance)*(fRmax + 0.5*kRadTolerance) ;
  tolIRMax2 = (fRmax - kRadTolerance*0.5)*(fRmax - kRadTolerance*0.5) ;

// Intersection with Rmax (possible return) and Rmin (must also check phi)

  G4int    i, j, num ;
  G4double Rtor2 = fRtor*fRtor, Rmax2 = fRmax*fRmax, Rmin2 = fRmin*fRmin ;
  G4double rho2  = p.x()*p.x()+p.y()*p.y();
  G4double rho   = sqrt(rho2) ;
  G4double pt2   = fabs(rho2+p.z()*p.z() +Rtor2 - 2*fRtor*rho) ;
  //   G4double pt = sqrt(pt2) ;
  G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
  G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;
  G4double vDotNmax = pDotV - fRtor*(v.x()*p.x() + v.y()*p.y())/rho ;

// Inside outer radius :
// check not inside, and heading through tubs (-> 0 to in)

  if( pt2 <= tolORMax2 && pt2 >= tolIRMin2 && vDotNmax < 0 )
  {
    if (seg)
    {
      inum   = p.x()*cosCPhi + p.y()*sinCPhi ;
      cosPsi = inum/rho ;

      if (cosPsi>=cosHDPhiIT) {
#if DEBUGTORUS
	G4cout << "G4Torus::DistanceToIn    (cosPsi>=cosHDPhiIT) "<< __LINE__ << G4endl << G4endl;
#endif	
	return snxt = 0 ;
      }
      
    }
    else {
#if DEBUGTORUS
      G4cout << "G4Torus::DistanceToIn    (seg) "<< __LINE__ << G4endl << G4endl;
#endif	
      return snxt = 0 ;
    }
  }
  else         // intersection with Rmax torus
  {    
    c[4] = 1.0 ;
    c[3] = 4*pDotV ;
    c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmax2 + 2*Rtor2*v.z()*v.z()) ;

    c[1] = 4*(pDotV*(pRad2 - Rtor2 - Rmax2) + 2*Rtor2*p.z()*v.z()) ;

    c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmax2) 
              + 4*Rtor2*p.z()*p.z() + (Rtor2-Rmax2)*(Rtor2-Rmax2) ;
   
    // num = SolveBiQuadratic(c,s) ;

    /* Numerical root research */
    s[0] = SolveNumeric(p, v, true);
    num = 1; // There is only one root: the correct one 
	
#if DEBUGTORUS
    G4cout << "G4Torus::DistanceToIn (" << __LINE__ << ") SolveNumeric : "
           << s[0] << G4endl;
#endif

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
	     xi     = p.x() + s[i]*v.x() ;
	     yi     = p.y() + s[i]*v.y() ;
 	     rhoi2  = xi*xi + yi*yi ;
             inum   = xi*cosCPhi + yi*sinCPhi ;
	     cosPsi = inum/sqrt(rhoi2) ;

	     if (cosPsi >= cosHDPhiIT)
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

    if( pt2 >= tolORMin2 && pt2 <= tolIRMax2 && vDotNmax > 0 )
    {
      if (seg)
      {
        inum   = p.x()*cosCPhi + p.y()*sinCPhi;
	cosPsi = inum/rho ;

	if (cosPsi>=cosHDPhiIT) {
#if DEBUGTORUS
	G4cout << "G4Torus::DistanceToIn    (cosPsi>=cosHDPhiIT) "<< __LINE__ << G4endl << G4endl;
#endif	
	  return snxt = 0 ;
	}
      }
      else {
#if DEBUGTORUS
	G4cout << "G4Torus::DistanceToIn     (seg) "<< __LINE__ << G4endl << G4endl;
#endif	
	return snxt = 0 ;
      }
    }
    else              // intersection with Rmin torus
    {               
      c[4] = 1.0 ;
      c[3] = 4*pDotV ;
      c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmin2 + 2*Rtor2*v.z()*v.z()) ;

      c[1] = 4*(pDotV*(pRad2-Rtor2-Rmin2) + 2*Rtor2*p.z()*v.z()) ;

      c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmin2) 
                    + 4*Rtor2*p.z()*p.z() + (Rtor2-Rmin2)*(Rtor2-Rmin2) ;
   
      // num = SolveBiQuadratic(c,s) ;

      /* Numerical root research */
      // s[0] = s[0]; // We already take care of Rmin in SolveNumeric !
      num = 1;

      if(num)
      {
        for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots   P?!
        {
     	  if(s[i] < kRadTolerance*0.5)
          {
	    for(j=i+1;j<num;j++) s[j-1] = s[j] ;
	    i-- ;
	    num-- ;
	  }
        }
	if(num)
	{
	  for(i = 0 ; i < num ; i++ )
	  {
	    if (seg)    // intersection point must have proper Phi
   	    {
	      xi     = p.x() + s[i]*v.x() ;
	      yi     = p.y() + s[i]*v.y() ;
 	      rhoi2  = xi*xi + yi*yi ;
              inum   = xi*cosCPhi + yi*sinCPhi ;
	      cosPsi = inum/sqrt(rhoi2) ;

	      if ( cosPsi >= cosHDPhiIT && s[i] < snxt )
	      {
	        snxt = s[i] ;
		break ;
	      }
	    }
	    else if(s[i] < snxt)
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

  if (seg)
  {                                      
    sinSPhi = sin(fSPhi) ; // First phi surface (`S'tarting phi)
    cosSPhi = cos(fSPhi) ;
    Comp    = v.x()*sinSPhi - v.y()*cosSPhi ;  // Compnent in outwards normal dirn
	                	
    if (Comp < 0 )
    {
      Dist = (p.y()*cosSPhi - p.x()*sinSPhi) ;

      if (Dist < kCarTolerance*0.5)
      {
        sphi = Dist/Comp ;

	if (sphi < snxt)
	{
	  if ( sphi < 0 ) sphi = 0 ;

          xi    = p.x() + sphi*v.x() ;
	  yi    = p.y() + sphi*v.y() ;
	  zi    = p.z() + sphi*v.z() ;
 	  rhoi2 = xi*xi + yi*yi ;
          it2   = fabs(rhoi2 + zi*zi + Rtor2 - 2*fRtor*sqrt(rhoi2)) ;

 	  if ( it2 >= tolORMin2 && it2 <= tolORMax2 )
 	  {
//  r intersection is good - check intersecting with correct half-plane

	    if ((yi*cosCPhi-xi*sinCPhi)<=0)	snxt=sphi;
 	  }    
	}
      }
    }
    ePhi=fSPhi+fDPhi;    // Second phi surface (`E'nding phi)
    sinEPhi=sin(ePhi);
    cosEPhi=cos(ePhi);
    Comp=-(v.x()*sinEPhi-v.y()*cosEPhi);
				
    if ( Comp < 0 )   // Component in outwards normal dirn
    {
      Dist = -(p.y()*cosEPhi - p.x()*sinEPhi) ;

      if (Dist < kCarTolerance*0.5 )
      {
        sphi = Dist/Comp ;

	if (sphi < snxt )
	{
	  if (sphi < 0 ) sphi = 0 ;
			 
          xi    = p.x() + sphi*v.x() ;
	  yi    = p.y() + sphi*v.y() ;
	  zi    = p.z() + sphi*v.z() ;
 	  rhoi2 = xi*xi + yi*yi ;
          it2   = fabs(rhoi2 + zi*zi + Rtor2 - 2*fRtor*sqrt(rhoi2)) ;

 	  if (it2 >= tolORMin2 && it2 <= tolORMax2)
 	  {
// z and r intersections good - check intersecting with correct half-plane

	    if ((yi*cosCPhi-xi*sinCPhi)>=0)	snxt=sphi;
 	  }    
	}
      }
    }
  }
  if(snxt < 0.5*kCarTolerance) snxt = 0.0 ;          

#if DEBUGTORUS
  G4cout << "G4Torus::DistanceToIn    Final Value is " << snxt << G4endl << G4endl;
#endif
  
  
  return snxt ;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Torus::DistanceToIn(const G4ThreeVector& p) const
{
  G4double safe, safe1, safe2 ;
  G4double phiC, cosPhiC, sinPhiC, safePhi, ePhi, cosPsi ;
  G4double rho2, rho, pt2, pt ;
    
#if DEBUGTORUS
	G4cout << G4endl ;
#endif

  rho2 = p.x()*p.x() + p.y()*p.y() ;
  rho  = sqrt(rho2) ;
  pt2  = fabs(rho2 + p.z()*p.z() + fRtor*fRtor - 2*fRtor*rho) ;
  pt   = sqrt(pt2) ;

  safe1 = fRmin - pt ;
  safe2 = pt - fRmax ;

  if (safe1 > safe2) safe = safe1;
  else               safe = safe2;

  if ( fDPhi < 2.0*M_PI && rho )
  {
    phiC    = fSPhi + fDPhi*0.5 ;
    cosPhiC = cos(phiC) ;
    sinPhiC = sin(phiC) ;
    cosPsi  = (p.x()*cosPhiC + p.y()*sinPhiC)/rho ;

    if (cosPsi < cos(fDPhi*0.5) ) // Psi=angle from central phi to point
    {                             // Point lies outside phi range
      if ((p.y()*cosPhiC - p.x()*sinPhiC) <= 0 )
      {
        safePhi = fabs(p.x()*sin(fSPhi) - p.y()*cos(fSPhi)) ;
      }
      else
      {
        ePhi    = fSPhi + fDPhi ;
	safePhi = fabs(p.x()*sin(ePhi) - p.y()*cos(ePhi)) ;
      }
      if (safePhi > safe) safe = safePhi ;
    }
  }
  if (safe < 0 ) safe = 0 ;
  return safe;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

/*
	Problem: if the ray exit the torus from the surface, the only solution is epsilon (~ 0). 
	Then this solution is eliminated by the loop (>= kRadTolerance) we have nothing.
	This results in 'invalid enum'
	solution: apply DistanceToIn instead DistanceToOut ?
 */

G4double G4Torus::DistanceToOut(const G4ThreeVector& p,
				const G4ThreeVector& v,
			        const G4bool calcNorm,
			        G4bool *validNorm,
				G4ThreeVector  *n    ) const
{
  ESide    side = kNull, sidephi = kNull ;
  G4double snxt = kInfinity, sphi, c[5], s[4] ;

  G4double sinSPhi, cosSPhi, ePhi, sinEPhi, cosEPhi;// Vars for phi intersection
  G4double cPhi, sinCPhi, cosCPhi ;

  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, zi, vphi ;


  // Radial Intersections Defenitions & General Precals
    
  // Define roots  Si (generally real >=0) for intersection with
  // torus (Ri = fRmax or fRmin) of ray p +S*v . General equation is :
  // c[4]*S^4 + c[3]*S^3 +c[2]*S^2 + c[1]*S + c[0] = 0 .
   
#if DEBUGTORUS
  G4cout << G4endl ;
#endif

  G4int    i,j,num ;
  G4double Rtor2 = fRtor*fRtor, Rmax2 = fRmax*fRmax, Rmin2 = fRmin*fRmin ;
  G4double rho2  = p.x()*p.x()+p.y()*p.y();
  G4double rho   = sqrt(rho2) ;
  G4double pt2   = fabs(rho2 + p.z()*p.z() + Rtor2 - 2*fRtor*rho) ;
  G4double pt    = sqrt(pt2) ;
  G4double pDotV = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;
  G4double pRad2 = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;
   
  G4double tolRMax = fRmax - kRadTolerance*0.5 ;
   
  G4double vDotNmax   = pDotV - fRtor*(v.x()*p.x() + v.y()*p.y())/rho ;
  G4double pDotxyNmax = (1 - fRtor/rho) ;

#if DEBUGTORUS
  G4cout << "G4Torus::DistanceToOut " << p << ", " << v << G4endl ;
#endif

  if( pt2 > tolRMax*tolRMax && vDotNmax >= 0 )
    {
      // On tolerant boundary & heading outwards (or perpendicular to) outer
      // radial surface -> leaving immediately with *n for really convex part only

      if (calcNorm && pDotxyNmax >= -kRadTolerance) 
	{
	  *n = G4ThreeVector( p.x()*(1 - fRtor/rho)/pt,
			      p.y()*(1 - fRtor/rho)/pt,
			      p.z()/pt                  ) ;
	  *validNorm = true ;
	}
#if DEBUGTORUS
      G4cout << "G4Torus::DistanceToOut    Leaving by Rmax immediately" << G4endl ;
#endif      
      return snxt = 0 ; // Leaving by Rmax immediately
    }
  else // intersection with Rmax torus
    {     
      c[4] = 1.0 ;
      c[3] = 4*pDotV ;
      c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmax2 + 2*Rtor2*v.z()*v.z()) ;

      c[1] = 4*(pDotV*(pRad2-Rtor2-Rmax2) + 2*Rtor2*p.z()*v.z()) ;

      c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmax2) 
	+ 4*Rtor2*p.z()*p.z() + (Rtor2-Rmax2)*(Rtor2-Rmax2) ;
   
				//num = SolveBiQuadratic(c,s) ;
				/* Numerical root research */
      s[0] = SolveNumeric( p, v, false);
      num = 1; // There is only one root.

#if DEBUGTORUS
      G4cout << "G4Torus::DistanceToOut (" << __LINE__ << ") SolveNumeric : " << s[0] << G4endl ;
#endif

      if(num)
	{
	  for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots
	    {
	      if( s[i] < kRadTolerance*0.5 )
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

      if (fRmin) // Possible Rmin intersection
	{
	  G4double tolRMin = fRmin + kRadTolerance*0.5 ;

	  // Leaving via Rmin
	  // NOTE: SHould use rho-rmin>kRadTolerance*0.5 - avoid sqrt for efficiency
	  if (pt2 < tolRMin*tolRMin && vDotNmax < 0 )
	    {
	      if (calcNorm)  *validNorm = false ;      // Concave surface of the torus
#if DEBUGTORUS
	      G4cout << "G4Torus::DistanceToOut    Leaving by Rmin immediately" << G4endl ;
#endif      
	      return          snxt      = 0 ;          // Leaving by Rmin immediately
	    }
	  else  // intersection with Rmin torus
	    {                
	      c[4] = 1.0 ;
	      c[3] = 4*pDotV ;
	      c[2] = 2*(pRad2 + 2*pDotV*pDotV - Rtor2 - Rmin2 + 2*Rtor2*v.z()*v.z()) ;

	      c[1] = 4*(pDotV*(pRad2-Rtor2-Rmin2) + 2*Rtor2*p.z()*v.z()) ;

	      c[0] = pRad2*pRad2 - 2*pRad2*(Rtor2+Rmin2) 
		+ 4*Rtor2*p.z()*p.z() + (Rtor2-Rmin2)*(Rtor2-Rmin2) ;
   
	      // num = SolveBiQuadratic(c,s) ;
	      // s[0] = s[0]; // We already take care of Rmin in SolveNumeric 
	      num = 1;
   
	      if(num)
		{
		  for(i=0;i<num;i++)   // leave only >=kRadTolerance/2 roots
		    {
		      if(s[i] < kRadTolerance*0.5)
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
  if (fDPhi < 2.0*M_PI)  // Phi Intersections
    {
      sinSPhi = sin(fSPhi) ;
      cosSPhi = cos(fSPhi) ;
      ePhi    = fSPhi + fDPhi ;
      sinEPhi = sin(ePhi) ;
      cosEPhi = cos(ePhi) ;
      cPhi    = fSPhi + fDPhi*0.5 ;
      sinCPhi = sin(cPhi) ;
      cosCPhi = cos(cPhi) ;


      if ( p.x() || p.y() ) // Check if on z axis (rho not needed later)
	{
	  pDistS = p.x()*sinSPhi - p.y()*cosSPhi ; // pDist -ve when inside
	  pDistE = -p.x()*sinEPhi + p.y()*cosEPhi ;

				// Comp -ve when in direction of outwards normal

	  compS   = -sinSPhi*v.x() + cosSPhi*v.y() ;
	  compE   = sinEPhi*v.x() - cosEPhi*v.y() ;
	  sidephi = kNull ;

	  if (pDistS <= 0 && pDistE <= 0 )
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

  /** Note: by numerical computation we know where the ray hits the torus **/
  /** So I propose to return the side where the ray hits **/
  if (calcNorm)
    {
      switch(side)
	{
	case kRMax:                     // n is unit vector 
#if DEBUGTORUS
	  G4cout << "G4Torus::DistanceToOut    Side is RMax" << G4endl ;
#endif
	  xi    = p.x() + snxt*v.x() ;
	  yi    =p.y() + snxt*v.y() ;
	  zi    = p.z() + snxt*v.z() ;
	  rhoi2 = xi*xi + yi*yi ;
	  rhoi  = sqrt(rhoi2) ;
	  it2   = fabs(rhoi2 + zi*zi + fRtor*fRtor - 2*fRtor*rhoi) ;
	  it    = sqrt(it2) ;
	  iDotxyNmax = (1-fRtor/rhoi) ;

	  if(iDotxyNmax >= -kRadTolerance) // really convex part of Rmax
	    {                       
	      *n = G4ThreeVector( xi*(1-fRtor/rhoi)/it,
				  yi*(1-fRtor/rhoi)/it,
				  zi/it                 ) ;
	      *validNorm = true ;
	    }
	  else *validNorm = false ; // concave-convex part of Rmax
	  break ;

	case kRMin:
#if DEBUGTORUS
	  G4cout << "G4Torus::DistanceToOut    Side is RMin" << G4endl ;
#endif
	  *validNorm = false ;	// Rmin is concave or concave-convex
	  break;

	case kSPhi:
#if DEBUGTORUS
	  G4cout << "G4Torus::DistanceToOut    Side is SPhi" << G4endl ;
#endif
	  if (fDPhi <= M_PI )
	    {
	      *n=G4ThreeVector(sin(fSPhi),-cos(fSPhi),0);
	      *validNorm=true;
	    }
	  else *validNorm = false ;
	  break ;

	case kEPhi:
#if DEBUGTORUS
	  G4cout << "G4Torus::DistanceToOut    Side is EPhi" << G4endl ;
#endif
	  if (fDPhi <= M_PI)
	    {
	      *n=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
	      *validNorm=true;
	    }
	  else *validNorm = false ;
	  break;

	default:

	  /* It seems we go here from time to time ..
	  G4cout << "Side is " << side << G4endl ;
	  G4cout << "Valid ESide are :" << kNull << " " << kRMin << " " << kRMax 
		 << " " << kSPhi << " " << kEPhi << G4endl;
	  */
	  G4cout.precision(16);
	  G4cout << G4endl;
	  G4cout << "Torus parameters:" << G4endl << G4endl;
	  G4cout << "fRmin = "   << fRmin/mm << " mm" << G4endl;
	  G4cout << "fRmax = "   << fRmax/mm << " mm" << G4endl;
	  G4cout << "fRtor = "   << fRtor/mm << " mm" << G4endl;
	  G4cout << "fSPhi = "   << fSPhi/degree << " degree" << G4endl;
	  G4cout << "fDPhi = "   << fDPhi/degree << " degree" << G4endl;
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
	
	  G4Exception("Invalid enum in G4Torus::DistanceToOut");
	  break;
	}
    }

#if DEBUGTORUS
  G4cout << "G4Torus::DistanceToOut    Final Value is " << snxt << G4endl << G4endl;
#endif

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
  rho2 = p.x()*p.x() + p.y()*p.y() ;
  rho  = sqrt(rho2) ;
  pt2  = fabs(rho2 + p.z()*p.z() + fRtor*fRtor - 2*fRtor*rho) ;
  pt   = sqrt(pt2) ;

#if DEBUGTORUS
	G4cout << G4endl ;
#endif

  if (fRmin)
  {
    safeR1 = pt - fRmin ;
    safeR2 = fRmax - pt ;

    if (safeR1 < safeR2) safe = safeR1 ;
    else                 safe = safeR2 ;
  }
  else safe = fRmax - pt ;    

// Check if phi divided, Calc distances closest phi plane

  if (fDPhi<2.0*M_PI) // Above/below central phi of Torus?
  {
    phiC    = fSPhi + fDPhi*0.5 ;
    cosPhiC = cos(phiC) ;
    sinPhiC = sin(phiC) ;

    if ((p.y()*cosPhiC-p.x()*sinPhiC)<=0)
    {
      safePhi = -(p.x()*sin(fSPhi) - p.y()*cos(fSPhi)) ;
    }
    else
    {
      ePhi    = fSPhi + fDPhi ;
      safePhi = (p.x()*sin(ePhi) - p.y()*cos(ePhi)) ;
    }
    if (safePhi < safe) safe = safePhi ;
  }
  if (safe < 0) safe = 0 ;
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
   G4Torus::CreateRotatedVertices(const G4AffineTransform& pTransform,
				  G4int& noPolygonVertices) const
{
  G4ThreeVectorList *vertices;
  G4ThreeVector vertex0,vertex1,vertex2,vertex3;
  G4double meshAngle,meshRMax,crossAngle,cosCrossAngle,sinCrossAngle,sAngle;
  G4double rMaxX,rMaxY,rMinX,rMinY;
  G4int crossSection,noCrossSections;

// Compute no of cross-sections necessary to mesh tube

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
  meshRMax  = (fRtor + fRmax)/cos(meshAngle*0.5) ;

// If complete in phi, set start angle such that mesh will be at fRmax
// on the x axis. Will give better extent calculations when not rotated.

  if ( fDPhi == M_PI*2.0 && fSPhi == 0 )
  {
    sAngle = -meshAngle*0.5 ;
  }
  else
  {
    sAngle = fSPhi ;
  }
  vertices = new G4ThreeVectorList(noCrossSections*4) ;
  
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

void G4Torus::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddThis (*this);
}

G4Polyhedron* G4Torus::CreatePolyhedron () const 
{
  return new G4PolyhedronTorus (fRmin, fRmax, fRtor, fSPhi, fSPhi + fDPhi);
}

G4NURBS* G4Torus::CreateNURBS () const 
{
  G4NURBS* pNURBS;
  if (fRmin != 0) 
  {
    if (fDPhi >= 2.0 * M_PI) 
    {
      pNURBS = new G4NURBStube (fRmin, fRmax, fRtor);
    }
    else 
    {
      pNURBS = new G4NURBStubesector (fRmin, fRmax, fRtor, fSPhi, fSPhi + fDPhi);
    }
  }
  else 
  {
    if (fDPhi >= 2.0 * M_PI) 
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


// 
// E.Medernach 
//

/** Important : the precision could be tuned by TORUSPRECISION **/

#define TORUSPRECISION 1.0  // or whatever you want for precision
                            // (it is TorusEquation related)
#define HOLEBVM 0
#define NBPOINT 6

/*
  Torus implementation with Newton Method and Bounding volume
 */

/*
  For speed issue,  we lose time *only* when intersecting the BVM
  and SafeNewton when it is called.
 */

G4double G4Torus::SolveNumeric(const G4ThreeVector& p,
                               const G4ThreeVector& v,
			       G4bool IsDistanceToIn) const
{
  /* This methods is a front-end to the numerical computation of roots */
  /* In fact this computation take care only of a perfect Torus */
  /* So here we add Phi section/Tolerance/Rinterior */

  /*** SolveNumeric ***/

  /** Conditions **/
  /** - if original point inside check for interior torus before **/
  /** - on surface it depends on the direction **/
  /** - the intersection point must be between fSPhi and fSPhi+fDPhi **/
  /** - If on the surface it depends on DistanceToOut or DistanceToIn : 
      a ray from the surface to In called with DistanceToIn return 0.0 
      and with DistanceToOut return the second intersection point **/

  G4double lambda = 0;
  G4double Value = TorusEquation(p.x(),p.y(),p.z(),GetRtor(),GetRmax());
  EInside inside ;


#if DEBUGTORUS
  G4cout << "G4Torus::SolveNumeric  " << p << ", " << v << G4endl ;
  G4cout << "G4Torus::SolveNumeric  Value = " << Value << G4endl;
#endif

  
  if (Value < -TORUSPRECISION) {
    inside = kInside ;
  } else {
    if (Value > TORUSPRECISION) {
      inside = kOutside;
    } else {
      inside = kSurface;
    }
  }
		

  switch (inside) {
  case kInside:
#if DEBUGTORUS
    G4cout << "G4Torus::SolveNumeric    Point is Inside Rmax Torus "
           << " Rtor = " << GetRtor()
	   << " Rmax = " << GetRmax() << G4endl ;
#endif
    if (fabs(GetRmin()) > POLEPSILON) {
#if DEBUGTORUS
      G4cout << "G4Torus::SolveNumeric    Testing interior torus .." << G4endl ;
#endif
      lambda = DistanceToTorus(p.x(),p.y(),p.z(),v.x(),v.y(),v.z(),
                               GetRtor(),GetRmin()); //Interior torus

#if DEBUGTORUS
      G4cout << "G4Torus::SolveNumeric    lambda to interior torus ="
             << lambda << G4endl ;
      G4cout << "G4Torus::SolveNumeric    Tolerance is "
             << kCarTolerance << G4endl ;
#endif
      /** Now check if on surface from interior torus **/

      /* PROBLEM: This may be a problem of precision
                  if we are near kCarTolerance ... */
      if (fabs(lambda) < kCarTolerance) {
	G4double Lx,Ly,Lz;
	G4double scal;
#if DEBUGTORUS								
	G4cout << "G4Torus::SolveNumeric    In fact on the Surface of Rmin torus"
	       << G4endl ;
#endif

	/* Compute Surface point */
	Lx = p.x() + lambda*v.x();
	Ly = p.y() + lambda*v.y();
	Lz = p.z() + lambda*v.z();
	/* Scalar product */
	scal  = v.x()*TorusDerivativeX(Lx,Ly,Lz,GetRtor(),GetRmin());
	scal += v.y()*TorusDerivativeY(Lx,Ly,Lz,GetRtor(),GetRmin());
	scal += v.z()*TorusDerivativeZ(Lx,Ly,Lz,GetRtor(),GetRmin());
	/* if entering and if it is DistToIn it is 0.0,  */
	/* but in fact it is the opposite because it is the interior torus */
	/* beware that this could be DistanceToOut */
	if ((IsDistanceToIn == true) && (scal > 0.0)) {
#if DEBUGTORUS
	  G4cout << "G4Torus::SolveNumeric    Entering Surface from Rmin Torus Gradient: "
	         << scal << G4endl ;
#endif
	  /* DistanceToIn return 0.0 */
	  lambda = 0.0;
	} else {
#if DEBUGTORUS
	  G4cout << "G4Torus::SolveNumeric    Exiting Surface (Recalculating) or DistanceToOut from surface"
	         << G4endl ;
	  G4cout << "G4Torus::SolveNumeric    Recursive call lambda..."
	         << lambda << G4endl << G4endl;
#endif
	  /* else it is not necessary infinity !!
	     (we could reach the opposite side..) */
	  /* To reach the opposite side we remark that from Surface
	     the sphere of radius min((Rmax - Rmin)/2, Rmin) does not
	     hit 2 surface of the torus so it is safe to do that way */
	  
	  if ((GetRmax() - GetRmin())/2.0 < GetRmin()) {
	    lambda = SolveNumeric(p+((GetRmax() - GetRmin())/2.0)*v,v,IsDistanceToIn) + (GetRmax() - GetRmin())/2.0;
	  } else {
	    lambda = SolveNumeric(p+GetRmin()*v,v,IsDistanceToIn) + GetRmin();
	  }

#if DEBUGTORUS
	  G4cout << "G4Torus::SolveNumeric    --> Recursive call: lambda = "
	         << lambda << G4endl;
#endif
	}
      } else {
	/* PROBLEM : could be better done ? */
	
	G4double lambdaToRmax = DistanceToTorus(p.x(),p.y(),p.z(),
	                                        v.x(),v.y(),v.z(),
						GetRtor(),GetRmax());
	if (lambda >= lambdaToRmax) {
#if DEBUGTORUS
	  G4cout << "G4Torus::SolveNumeric    Point does not hit the Rmin torus from here" << G4endl;
#endif
	  lambda = lambdaToRmax; 
	} else {
#if DEBUGTORUS							
	  G4cout << "G4Torus::SolveNumeric    We hit the Rmin torus with "
	         << lambda << G4endl;
	  G4cout << "G4Torus::SolveNumeric    Note that this could be small and not in Tolerance resulting in wrong result " 
		 << G4endl ;
#endif
	}
      }
    } else {
      /* It is a whole torus */
      
      lambda = DistanceToTorus(p.x(),p.y(),p.z(),
                               v.x(),v.y(),v.z(),
			       GetRtor(),GetRmax()); 
    }
    break;
  case kSurface:
    {
      G4double Lx,Ly,Lz;
      G4double scal;

#if DEBUGTORUS
      G4cout << "G4Torus::SolveNumeric    Point is on the Rmax Surface"
             << G4endl ;
#endif
      /* It is possible with Phi that this is not the correct point */
      lambda = DistanceToTorus(p.x(),p.y(),p.z(),
                               v.x(),v.y(),v.z(),
			       GetRtor(),GetRmax()); 
      /* Compute Surface point */
      Lx = p.x() + lambda*v.x();
      Ly = p.y() + lambda*v.y();
      Lz = p.z() + lambda*v.z();
      /* Scalar product */
      scal  = v.x()*TorusDerivativeX(Lx,Ly,Lz,GetRtor(),GetRmax());
      scal += v.y()*TorusDerivativeY(Lx,Ly,Lz,GetRtor(),GetRmax());
      scal += v.z()*TorusDerivativeZ(Lx,Ly,Lz,GetRtor(),GetRmax());
      
      /* if entering it is < 0.0 */
      if ((IsDistanceToIn) && (scal < 0.0)) {
#if DEBUGTORUS
	G4cout << "G4Torus::SolveNumeric    Point is Entering Surface "
	       << scal << G4endl ;
#endif
	lambda = 0.0;
      } else {
#if DEBUGTORUS
	G4cout << "G4Torus::SolveNumeric    Point is Exiting Surface or DistanceToOut "
	       << scal << G4endl ;
	G4cout << "Recursive call ..." << G4endl << G4endl ;
#endif
	/* To reach the opposite side we remark that from Surface the sphere of radius (Rmax - Rmin)/2 
	   does not hit 2 surface of the torus so it is safe to do that way */
	//lambda = SolveNumeric(p+(lambda + kCarTolerance)*v,v,IsDistanceToIn);
	lambda = SolveNumeric(p+((GetRmax() - GetRmin())/2.0)*v,
	                      v, IsDistanceToIn)
		 + (GetRmax() - GetRmin())/2.0;
#if DEBUGTORUS
	G4cout << "Recursive call ...END" << G4endl ;
#endif
      }
    }	
    break;
  case kOutside:
#if DEBUGTORUS
    G4cout << "G4Torus::SolveNumeric    Point is Outside the Rmax torus"
           << G4endl ;
#endif
	       
    lambda = DistanceToTorus(p.x(),p.y(),p.z(),
                             v.x(),v.y(),v.z(),
			     GetRtor(),GetRmax()); 
    break;
  }

  if (lambda == kInfinity) return lambda;

#if DEBUGTORUS
  G4cout << "G4Torus::SolveNumeric    Intersection found. Now checking Phi angles"
         << G4endl ;
#endif
  
  /** Ok we have a lambda that is correct without Phi **/
  /** Now check Phi .. **/

  /* Eliminate the case of point (0,0,0) */
  //  if ((p.x()*p.x() + p.y()*p.y() + p.z()*p.z()) > POLEPSILON)
  if (((p.x()+ lambda*v.x())*(p.x()+ lambda*v.x()) +
       (p.y()+ lambda*v.y())*(p.y()+ lambda*v.y()) +
       (p.z()+ lambda*v.z())*(p.z()+ lambda*v.z())) > POLEPSILON)
  {
    G4double theta ;

    theta = atan2(p.y() + lambda*v.y(),p.x() + lambda*v.x());
#if DEBUGTORUS
    G4cout << "G4Torus::SolveNumeric    theta = " << theta << G4endl;
#endif 

    if (theta < 0) theta += 2*M_PI;
    
    
    /*** We have to verify if this root is inside the region between fSPhi and fSPhi + fDPhi ***/
#if DEBUGTORUS
    G4cout << "G4Torus::SolveNumeric    theta = " << theta
           << " Phi = " << fSPhi 
	   << " Phi + dPhi = " << fSPhi + fDPhi
	   << " kAngTolerance = " << kAngTolerance << G4endl ;
    G4cout << " theta - Phi = " << theta - fSPhi << G4endl ;
    
#endif 
    
    if ((theta - fSPhi >= - kAngTolerance*0.5) &&
        (theta - (fSPhi + fDPhi) <=  kAngTolerance*0.5)) {
      /*** If this is the case we return this solution ***/
#if DEBUGTORUS
      G4cout << "G4Torus::SolveNumeric    Correct Phi section" << G4endl ;
#endif
      return lambda;
    } else {
      /*** Else we compute the intersection with the 2 half-plane [fSPhi]
           and [fSPhi + fDPhi] ***/

      G4double IntersectPlanar ;

      
      IntersectPlanar = - (p.y() - p.x()*tan(fSPhi))/(v.y() - v.x()*tan(fSPhi));
#if DEBUGTORUS
      G4cout << "G4Torus::SolveNumeric    IntersectPlanar = "
             << IntersectPlanar << G4endl ;
#endif

      /** If this is below lambda we check for the other plane **/
      if (IntersectPlanar < lambda) { 
	IntersectPlanar = - (p.y() - p.x()*tan(fSPhi + fDPhi))
	                  / (v.y() - v.x()*tan(fSPhi + fDPhi)) ;
#if DEBUGTORUS
	G4cout << "G4Torus::SolveNumeric    IntersectPlanar (2) = "
	       << IntersectPlanar << G4endl ;
#endif
      }
      
      /* If we does not hit the two plan then we does not hit the torus .. */
      if (IntersectPlanar < lambda) {
#if DEBUGTORUS
	G4cout << "G4Torus::SolveNumeric    No intersection with planar Phi .."
	       << G4endl ;
#endif
	  return kInfinity;
      }
      
#if DEBUGTORUS
      G4cout << "G4Torus::SolveNumeric    Incorrect Phi section" << G4endl ;
      G4cout << "G4Torus::SolveNumeric    point : " << p << " direction : "
             << v << G4endl ;
      G4cout << "G4Torus::SolveNumeric    IntersectPlanar = "
             << IntersectPlanar << G4endl ;
#endif
      
      if ((TorusEquation(p.x() + IntersectPlanar*v.x(),
			 p.y() + IntersectPlanar*v.y(),
			 p.z() + IntersectPlanar*v.z(),
			 GetRtor(),GetRmax()) < 0)
	  &&  (TorusEquation(p.x() + IntersectPlanar*v.x(),
			     p.y() + IntersectPlanar*v.y(),
			     p.z() + IntersectPlanar*v.z(),
			     GetRtor(),GetRmin()) > 0)) {
	/*** if this point is inside torus Rmax and outside torus Rmin
	     then it is on the cut planar faces ***/
#if DEBUGTORUS
	G4cout << "G4Torus::SolveNumeric    Hit planar section" << G4endl ;
#endif
	return IntersectPlanar;
      } else {
	/*** else we continue from this new point (SolveNumeric) ***/
#if DEBUGTORUS
	G4cout << "G4Torus::SolveNumeric    Recursive Phi call with "
	       << IntersectPlanar << " .." << G4endl << G4endl;
#endif

	return IntersectPlanar + SolveNumeric(p+IntersectPlanar*v,
	                                      v,IsDistanceToIn);
      }
    }
  } else {
#if DEBUGTORUS
    G4cout << "G4Torus::SolveNumeric    Phi not checked because point is " << p + lambda*v << G4endl << G4endl;
#endif
  }
  

  return lambda;
}

void G4Torus::BVMIntersection(G4double x,G4double y,G4double z,
			      G4double dx,G4double dy,G4double dz,
			      G4double Rmax, G4double Rmin,
			      G4double *NewL,G4int *valid) const
{

  if (dz != 0) {
    G4double DistToZ ;
    /* z = + Rmin */
    NewL[0] = (Rmin - z)/dz ;
    /* z = - Rmin */
    NewL[1] = (-Rmin - z)/dz ;
    /* Test validity here (*** To be optimized ***) */
    if (NewL[0] < 0.0) valid[0] = 0;
    if (NewL[1] < 0.0) valid[1] = 0;
    DistToZ = (x+NewL[0]*dx)*(x+NewL[0]*dx) + (y+NewL[0]*dy)*(y+NewL[0]*dy);
    if (DistToZ	- (Rmax + Rmin)*(Rmax + Rmin) > 0)
      valid[0] = 0;
#if HOLEBVM 
    if (DistToZ	- (Rmax - Rmin)*(Rmax - Rmin) < 0)
      valid[0] = 0;
#endif
    DistToZ = (x+NewL[1]*dx)*(x+NewL[1]*dx) + (y+NewL[1]*dy)*(y+NewL[1]*dy);
    if (DistToZ	- (Rmax + Rmin)*(Rmax + Rmin) > 0)
      valid[1] = 0;
#if HOLEBVM
    if (DistToZ	- (Rmax - Rmin)*(Rmax - Rmin) < 0)
      valid[1] = 0;
#endif    
  } else {
    /* if dz == 0 we could know the exact solution */
    /* Well, this is true but we have not expected precision issue from sqrt .. */
    NewL[0] = -1.0;
    NewL[1] = -1.0;
    valid[0] = 0;
    valid[1] = 0;
  }

  /* x� + y� = (Rmax + Rmin)� */
  if ((dx != 0) || (dy != 0)) {
    G4double a,b,c,d;
    
    a = dx*dx + dy*dy ;
    b = 2*(x*dx + y*dy) ;
    c = x*x + y*y - (Rmax + Rmin)*(Rmax + Rmin) ;
    d = b*b - 4*a*c ;
    
    if (d < 0) {
      valid[2] = 0;
      valid[3] = 0;
      NewL[2] = -1.0;
      NewL[3] = -1.0;
    } else{
      d = sqrt(d) ;
      NewL[2] = (d - b)/(2*a);
      NewL[3] = (-d - b)/(2*a);
      if (NewL[2] < 0.0) valid[2] = 0;
      if (fabs(z + NewL[2]*dz) - Rmin > POLEPSILON) valid[2] = 0;
      if (NewL[3] < 0.0) valid[3] = 0;
      if (fabs(z + NewL[3]*dz) - Rmin > POLEPSILON) valid[3] = 0;
    }
  } else
    {
      /* only dz != 0 so we could know the exact solution */
      /* this depends only for the distance to Z axis */
      /* BUT big precision problem near the border.. */
      /* I like so much Newton to increase precision you know.. */

      NewL[2] = -1.0;
      NewL[3] = -1.0;
      valid[2] = 0;
      valid[3] = 0;
	
      /*** Try This to see precision issue with sqrt(~ 0)
	   G4double DistToZ ;
	   G4double result;
	   G4double guess;
	
	   DistToZ = sqrt(x*x + y*y) ;
	
	   if ((DistToZ < (Rmax - Rmin)) || (DistToZ > (Rmax + Rmin))) {
	   return -1.0 ;
	   }
	
	   result = sqrt((Rmin + Rmax - DistToZ)*(Rmin - Rmax + DistToZ));

	   if (dz < 0) {
	   if (z > result) {
	   return (result - z)/dz;
	   } else {
	   if (z > -result) {
	   return (-result - z)/dz;
	   } else 
	   return -1.0;
	   }
	   } else {
	   if (z < -result) {
	   return (z + result)/dz;
	   } else {
	   if (z < result) {
	   return (z - result)/dz;
	   } else 
	   return -1.0;
	   }
	   }
      */
    }
  

  /* x� + y� = (Rmax - Rmin)� */
#if HOLEBVM
  if ((dx != 0) || (dy != 0)) {
    G4double a,b,c,d;
    
    a = dx*dx + dy*dy ;
    b = 2*(x*dx + y*dy) ;
    c = x*x + y*y - (Rmax - Rmin)*(Rmax - Rmin) ;
    d = b*b - 4*a*c ;
    
    if (d < 0) {
      valid[4] = 0;
      valid[5] = 0;
      NewL[4] = -1.0;
      NewL[5] = -1.0;
    } else {
      d = sqrt(d) ;
      NewL[4] = (d - b)/(2*a);
      NewL[5] = (-d - b)/(2*a);
      if (NewL[4] < 0.0) valid[4] = 0;
      if (fabs(z + NewL[4]*dz) - Rmin > POLEPSILON) valid[4] = 0;
      if (NewL[5] < 0.0) valid[5] = 0;
      if (fabs(z + NewL[5]*dz) - Rmin > POLEPSILON) valid[5] = 0;
    }
  } else
#endif            
    {
      /* only dz != 0 so we could know the exact solution */
      /* OK but same as above .. */
      valid[4] = 0;
      valid[5] = 0;
      NewL[4] = -1.0;
      NewL[5] = -1.0;
    }
}

void G4Torus::SortIntervals (G4double *SortL, G4double *NewL,
                             G4int *valid, G4int *NbIntersection) const
{
  G4int i,j;
  G4double swap;
	
  (*NbIntersection) = 0;
  SortL[0] = -kInfinity;
	
  for (i=0;i<6;i++) {
    if (valid[i] != 0) {
      SortL[(*NbIntersection)] = NewL[i] ;
      for (j=(*NbIntersection);j>0;j--) {
	if (SortL[j] < SortL[j-1]) {
	  swap = SortL[j-1] ;
	  SortL[j-1] = SortL[j];
	  SortL[j] = swap;
	}
      }
		
      (*NbIntersection) ++;
    }
  }
  /* Delete double values */
  /* When the ray hits a corner we have a double value */
  for (i=0;i<(*NbIntersection)-1;i++) {
    if (SortL[i+1] - SortL[i] < POLEPSILON) {
      if (((*NbIntersection) & (1)) == 1) {
	/* If the NbIntersection is odd then we keep one value */
	for (j=i+1;j<(*NbIntersection);j++) {
	  SortL[j-1] = SortL[j] ;
	}
	(*NbIntersection) --;
      } else {
	/* If it is even we delete the 2 values */
	for (j=i+2;j<(*NbIntersection);j++) {
	  SortL[j-2] = SortL[j] ;
	}
	(*NbIntersection) -= 2;
      }
    }
  }
}



/** Now the interesting part .. **/


G4double G4Torus::DistanceToTorus (G4double x,G4double y,G4double z,
				   G4double dx,G4double dy,G4double dz,
				   G4double Rmax,G4double Rmin) const
{
  G4double Lmin=0.;
  G4double Lmax=0.;
  G4double guess;
  G4double SortL[4];
   
  G4int NbIntersection = 0;

  G4double NewL[NBPOINT];
  G4int valid[] = {1,1,1,1,1,1} ;
  G4int j;

  j = 0;


  /*** Compute Intervals  from Bounding Volume ***/

  BVMIntersection(x,y,z,dx,dy,dz,Rmax,Rmin,NewL,valid);

  /*
    We could compute intervals value 
    Sort all valid NewL to SortL.
    There must be 4 values at max and 
    odd one if point is inside
  */

  SortIntervals(SortL,NewL,valid,&NbIntersection);
  
  {
    /*** Length check (Torus specific) ***/
    G4double LengthMin = 0.82842712*Rmin;
				
    switch(NbIntersection) {
    case 1:
      if (SortL[0] < POLEPSILON) {
	if (fabs(TorusEquation(x,y,z,Rmax,Rmin)) < TORUSPRECISION) {
	  return 0.0;
	} else {
	  return kInfinity;
	}
      }
      break;
    case 2:
      if ((SortL[1] - SortL[0]) < LengthMin) NbIntersection = 0;
      break;
    case 3:
      if (SortL[0] < POLEPSILON) {
	if (fabs(TorusEquation(x,y,z,Rmax,Rmin)) < TORUSPRECISION) {
	  return 0.0;
	} else {
	  NbIntersection --;
	  SortL[0] = SortL[1] ;
	  SortL[1] = SortL[2] ;
	  if ((SortL[1] - SortL[0]) < LengthMin) NbIntersection = 0;
	}
      } else {
	if ((SortL[2] - SortL[1]) < LengthMin) NbIntersection -= 2;
      }
      break;
    case 4:
      if ((SortL[1] - SortL[0]) < LengthMin) {
	NbIntersection -= 2;
	SortL[0] = SortL[2];
	SortL[1] = SortL[3];
	if ((SortL[1] - SortL[0]) < LengthMin) NbIntersection -= 2;	
      }
      break;
    }
  }
	
#if DEBUGTORUS
  {
    G4int i;
    G4cout.precision(16);
    G4cout << "G4Torus::DistanceToTorus    INTERVALS" << G4endl ;
    for (i=0;i<NbIntersection;i++) {
      G4cout << "G4Torus::DistanceToTorus    " << SortL[i] << G4endl ;
    }
  }
#endif

  switch (NbIntersection) {
  case 0:
    return kInfinity ;				
    break;
  case 1:
    Lmin = 0.0 ;
    Lmax  = SortL[0] ;
    break;
  case 2:
    Lmin = SortL[0] ;
    Lmax = SortL[1] ;
    break;
#if HOLEBVM
  case 3:
    Lmin = 0.0 ;
    Lmax = SortL[0] ;
    
    TorusEquationClass torus (Rmax,Rmin);
    torus.setPosition(x,y,z);
    torus.setDirection(dx,dy,dz);
  
    G4PolynomialSolver<TorusEquationClass,G4double(TorusEquationClass::*)(G4double)>
      PolySolver(&torus,
		 &TorusEquationClass::Function,
		 &TorusEquationClass::Derivative,
		 TORUSPRECISION) ;

    guess = PolySolver.solve(Lmin,Lmax);

    if ((guess >= (Lmin - POLEPSILON)) && (guess <= (Lmax + POLEPSILON))) {
      return guess ;
    } else {
      Lmin = SortL[1] ;
      Lmax = SortL[2] ;
    }
    break;
  case 4:
    Lmin = SortL[0] ;
    Lmax = SortL[1] ;

    TorusEquationClass torus (Rmax,Rmin);
    torus.setPosition(x,y,z);
    torus.setDirection(dx,dy,dz);
  
    G4PolynomialSolver<TorusEquationClass,G4double(TorusEquationClass::*)(G4double)>
      PolySolver(&torus,
		 &TorusEquationClass::Function,
		 &TorusEquationClass::Derivative,
		 TORUSPRECISION) ;

    guess = PolySolver.solve(Lmin,Lmax);

    if ((guess >= (Lmin - POLEPSILON)) && (guess <= (Lmax + POLEPSILON))) {
      return guess ;
    } else {
      Lmin = SortL[2] ;
      Lmax = SortL[3] ;
    }
    break;
#endif
    
  default:
    G4cerr << "G4Torus::DistanceToTorus    NbIntersection = " << NbIntersection << G4endl;    
    break;    
  }

  TorusEquationClass torus (Rmax,Rmin);
  torus.setPosition(x,y,z);
  torus.setDirection(dx,dy,dz);
  
  G4PolynomialSolver<TorusEquationClass,G4double(TorusEquationClass::*)(G4double)>
    PolySolver(&torus,
	       &TorusEquationClass::Function,
	       &TorusEquationClass::Derivative,
	       TORUSPRECISION) ;

  guess = PolySolver.solve(Lmin,Lmax);

  if ((guess >= (Lmin - POLEPSILON)) && (guess <= (Lmax + POLEPSILON))) {
#if DEBUGTORUS
    G4cout << "G4Torus::DistanceToTorus    distance = " << guess << G4endl ;    
#endif
    return guess ;
  } else {
#if DEBUGTORUS
    G4cout << "G4Torus::DistanceToTorus  :  kInfinity" << G4endl ;    
#endif
    return kInfinity;
  }
}
