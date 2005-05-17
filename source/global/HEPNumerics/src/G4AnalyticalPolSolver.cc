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
// $Id: G4AnalyticalPolSolver.cc,v 1.2 2005-05-17 14:49:21 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4AnalyticalPolSolver.hh"

////////////////////////////////////////////

G4AnalyticalPolSolver::G4AnalyticalPolSolver() {;}

////////////////////////////////////////////////////////

G4AnalyticalPolSolver::~G4AnalyticalPolSolver() {;}

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^2 + p[1] x+p[2]=0
//
// x = r[1][k] + i r[2][k];  k = 1, 2

G4int G4AnalyticalPolSolver::QuadRoots( G4double p[5], G4double r[3][5] )
{
  G4double b, c, d, tmp;

  b = -p[1]/p[0]/2;
  c =  p[2]/p[0];
  d =  b*b - c;

  if( d > 0 )
  {
    if( b > 0 ) b = r[1][2] = b + std::sqrt(d);
    else        b = r[1][2] = b - std::sqrt(d);

    r[1][1] = c/b; 
    r[2][1] = 0.; 
    r[2][2] = 0.;
  }
  else
  {
    tmp     =  std::sqrt(-d);
    r[2][1] =  tmp; 
    r[2][2] = -tmp;
    r[1][1] =  b; 
    r[1][2] =  b;
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^3 + p[1] x^2...+p[3]=0
// x=r[1][k] + i r[2][k]  k=1,...,3
// Assumes 0<arctan(x)<pi/2 for x>0

G4int G4AnalyticalPolSolver::CubicRoots( G4double p[5], G4double r[3][5] )
{
  G4double s,t,b,c,d;
  G4int k;

  if( p[0] != 1. )
  {
   for(k = 1; k < 4; k++ ) p[k] = p[k]/p[0]; 
                           p[0] = 1.;
  }
  s = p[1]/3.0; 
  t = s*p[1];
  b = 0.5*( s*( t/1.5 - p[2] ) + p[3] ); 
  t = ( t - p[2] )/3.0;
  c = t*t*t; 
  d = b*b - c;

  if( d >= 0. )
  {
    d = std::pow( (std::sqrt(d) + std::fabs(b) ), 1.0/3.0 );
    
    if( d != 0. )
    {
       if( b > 0. ) b = -d;
       else         b =  d;
                    c =  t/b;
    }
    d       =  std::sqrt(0.75)*(b - c); 
    r[2][2] =  d; 
    b       =  b + c;
    c       = -0.5*b-s;
    r[1][2] =  c;

    if( ( b > 0. &&  s <= 0. ) || ( b < 0. && s > 0. ) )
    {
       r[1][1] =  c; 
       r[2][1] = -d; 
       r[1][3] =  b - s;
       r[2][3] =  0;
    }
    else
    {
       r[1][1] =  b - s; 
       r[2][1] =  0.; 
       r[1][3] =  c;
       r[2][3] = -d;
    }
  }              // end of 2 equal or complex roots 
  else
  {
    if( b == 0. ) d =  std::atan(1.0)/1.5;
    else          d =  std::atan( std::sqrt(-d)/std::fabs(b) )/3.0;

    if( b < 0. )  b =  std::sqrt(t)*2.0;
    else          b = -2.0*std::sqrt(t);

    c =  std::cos(d)*b; 
    t = -std::sqrt(0.75)*std::sin(d)*b - 0.5*c;
    d = -t - c - s; 
    c =  c - s; 
    t =  t - s;

    if( std::fabs(c) > std::fabs(t) ) r[1][3] = c;
    else
    {
      r[1][3] = t; 
      t       = c;
    }
    if( std::fabs(d) > std::fabs(t) ) r[1][2] = d;
    else
    {
      r[1][2] = t; 
      t       = d;
    }
    r[1][1] = t;

    for(k = 1; k < 4; k++ ) r[2][k] = 0. ;
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
// x=r[1][k] + i r[2][k]  k=1,...,4

G4int G4AnalyticalPolSolver::BiquadRoots( G4double p[5], G4double r[3][5] )
{
  G4double a, b, c, d, e;
  G4int k, j;

  if(p[0] != 1.0)
  {
    for( k = 1; k < 5; k++) p[k] = p[k]/p[0];
                            p[0] = 1.;
  }
  e = 0.25*p[1];
  b = 2*e;
  c = b*b;
  d = 0.75*c;
  b = p[3] + b*( c - p[2] );
  a = p[2] - d;
  c = p[4] + e*( e*a - p[3] );
  a = a - d;

  p[1] = 0.5*a;
  p[2] = (p[1]*p[1]-c)*0.25;
  p[3] = b*b/(-64.0);
  if( p[3] < -1.e-6 )
  {
    CubicRoots(p,r);

    for( k = 1; k < 4; k++ )
    {
      if( r[2][k] == 0. && r[1][k] > 0 )
      {
        d = r[1][k]*4; 
        a = a + d;

        if     ( a >= 0. && b >= 0.) p[1] =  std::sqrt(d);
        else if( a <= 0. && b <= 0.) p[1] =  std::sqrt(d);
        else                         p[1] = -std::sqrt(d);

        b = 0.5*( a + b/p[1] );
        break ;
      }
    }
  }
  if( p[2] < 0. )
  {
    b    = std::sqrt(c); 
    d    = b + b - a;
    p[1] = 0.; 
    
    if( d > 0. ) p[1] = std::sqrt(d);
  }
  else
  {
    if( p[1] > 0.) b =  std::sqrt(p[2])*2.0 + p[1];
    else           b = -std::sqrt(p[2])*2.0 + p[1];

    if(b != 0.) p[1] = 0;
    else
    {
      for(k = 1; k < 5; k++ )
      {
          r[1][k] = -e;
          r[2][k] =  0;
      }
      return 0;
    }
  }

  p[2] = c/b; 
  QuadRoots(p,r);

  for( k = 1; k < 3; k++ )
  {
    for( j = 1; j < 3; j++ ) r[j][k+2] = r[j][k];
  }
  p[1] = -p[1]; 
  p[2] =  b; 
  QuadRoots(p,r);

  for( k =1 ; k < 5; k++ ) r[1][k] = r[1][k] - e;

  return 0;    //  return(0);
}

//
//
//////////////////////////////////////////////////////////////////////////////
