//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4ApproxPolySolver.cc
//
// Author:
//
//   Oliver Link (oliver.Link@cern.ch)
// --------------------------------------------------------------------

#include "G4ApproxPolySolver.hh" 

/////////////////////////////////////////////////////////////////////////
//
// Auxiliary method for solving (in real numbers) biquadratic equation
// Algorithm based on : Graphics Gems I by Jochen Schwartz

G4int G4ApproxPolySolver::SolveBiQuadratic( G4double c[], G4double s[]  ) const
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

  // y^4 + py^2 + r = 0 and z=y^2 so y = +-std::sqrt(z1) and y = +-std::sqrt(z2)
   
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
            s[2] = std::sqrt(s[1]) ;
            s[1] = s[0] ;
            s[0] = -s[2] ;
            num++ ;
          }
          else        // Four roots
          {
            s[2] = std::sqrt(s[0]) ;
            s[3] = std::sqrt(s[1]) ;
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
            s[0] = -std::sqrt(s[1]) ;
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
            s[1] = std::sqrt(s[0]) ;
            s[0] = -s[1] ;
            num +=1 ;
          }
        }
        else return num = 0 ;
      }
    }
    else return num ;
  }
  else if (r == 0)     // no absolute term: y(y^3 + py + q) = 0 
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
    else if (u > 0)  u = std::sqrt(u) ;
    else             return 0 ;

    if (v==0)        v = 0 ;
    else if (v > 0)  v = std::sqrt(v);
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

G4int G4ApproxPolySolver::SolveCubic( G4double c[], G4double s[] ) const
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
  //  x^3 +px + q = 0 

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
      G4double u = std::pow(-q,1./3.);
      s[ 0 ] = 2 * u;
      s[ 1 ] = - u;
      num = 2;
    }
  }
  else if (D < 0) // Casus irreducibilis: three real solutions
  {
    G4double phi = 1.0/3 * std::acos(-q / std::sqrt(-p3));
    G4double t = 2 * std::sqrt(-p);

    s[ 0 ] =   t * std::cos(phi);
    s[ 1 ] = - t * std::cos(phi + pi / 3);
    s[ 2 ] = - t * std::cos(phi - pi / 3);
    num = 3;
  }
  else // one real solution 
  {
    G4double sqrt_D = std::sqrt(D);
    G4double u = std::pow(sqrt_D - q,1./3.);
    G4double v = - std::pow(sqrt_D + q,1./3.);

    s[ 0 ] = u + v;
    num = 1;
  }

  // resubstitute 

  sub = 1.0/3 * A;

  for (i = 0; i < num; ++i)
    s[ i ] -= sub;

  return num;
}

////////////////////////////////////////////////////////////////////////////
//
//

G4int G4ApproxPolySolver::SolveBiQuadraticNew( G4double c[], G4double s[] ) const
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
      s[0] = std::sqrt( std::sqrt( -D ) );
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
      v[i] = std::fabs( s[i] ) ;
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

    w1 = std::sqrt( s[i] );
    w2 = std::sqrt( s[j] );
  } 
  else 
  {
    num = 2;
    w1 = w2 = std::sqrt( s[1] );
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

///////////////////////////////////////////////////////////////////////////
//
//

G4int G4ApproxPolySolver::SolveCubicNew( G4double c[], G4double s[],
                              G4double& cubic_discr ) const
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
    h2 = std::sqrt( cubic_discr );
    u = -h1+h2;
    v = -h1-h2;
    if( u < 0 ) u = -std::pow(-u,1./3.);
    else u = std::pow(u,1./3.);
    if( v < 0 ) v = -std::pow(-v,1./3.);
    else v = std::pow(v,1./3.);
    s[0] = u+v-sub;
    s[1] = -(u+v)/2.0-sub;
    s[2] = std::fabs(u-v)*std::sqrt(3.0)/2.0;
    if( std::fabs(u) <= eps || std::fabs(v) <= eps ) 
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
  else if( std::fabs(cubic_discr) <= delta ) 
  {
    cubic_discr = 0.;

    if( h1 < 0 ) u = std::pow(-h1,1./3.);
    else         u = -std::pow(h1,1./3.);

    s[0] =  u + u - sub ;
    s[1] = -u - sub ;
    s[2] = s[1] ;

    if( std::fabs(h1) <= eps ) 
    {
      y[0] = s[0];
      for( i=0; i<2; i++ ) 
      {
        h1 = (3.0*y[i]+2.*A)*y[i]+B;

        if( std::fabs(h1) > delta )
          y[i+1] = y[i]-(((y[i]+A)*y[i]+B)*y[i]+C)/h1;
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
    h3 =std::fabs(p/3.);
    h3 = std::sqrt(h3*h3*h3);
    h2 = std::acos(-h1/h3)/3.;
    h1 = std::pow(h3,1./3.);
    u = h1*std::cos(h2);
    v = std::sqrt(3.)*h1*std::sin(h2);
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

G4int G4ApproxPolySolver::SolveQuadratic( G4double c[], G4double s[] ) const
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
    G4double sqrt_D = std::sqrt(D);

    s[ 0 ] = - p - sqrt_D ;  // in ascending order !
    s[ 1 ] = - p + sqrt_D ;
    return 2;
  }
  return 0;
}
