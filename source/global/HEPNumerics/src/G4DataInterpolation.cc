// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DataInterpolation.cc,v 1.2 1999-08-24 12:00:12 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4DataInterpolation.hh"


// Constructor for initializing of fArgument, fFunction and fNumber data members

G4DataInterpolation::G4DataInterpolation( G4double pX[], 
	                                  G4double pY[], 
			                  G4int number    ):
  fSecondDerivative (0)
{
   G4int i ;
   fNumber = number ;
   fArgument = new G4double[fNumber] ;
   fFunction = new G4double[fNumber] ;
   for(i=0;i<fNumber;i++)
   {
      fArgument[i] = pX[i] ;
      fFunction[i] = pY[i] ;
   }
} 

// Constructor for cubic spline interpolation. It creates the array 
// fSecondDerivative[0,...fNumber-1] which is used in this interpolation by
// the function 


G4DataInterpolation::G4DataInterpolation( G4double pX[], 
	                                  G4double pY[], 
			                  G4int number,
				          G4double pFirstDerStart,
				          G4double pFirstDerFinish  ) 
{
   G4int i, k ;
   G4double p, qn, sig, un ;
   const G4double maxDerivative = 0.99e30 ;
   fNumber = number ;
   fArgument = new G4double[fNumber] ;
   fFunction = new G4double[fNumber] ;
   fSecondDerivative = new G4double[fNumber] ;
   G4double* u = new G4double[fNumber - 1] ;

   for(i=0;i<fNumber;i++)
   {
      fArgument[i] = pX[i] ;
      fFunction[i] = pY[i] ;
   }
   if(pFirstDerStart > maxDerivative)
   {
      fSecondDerivative[0] = 0.0 ;
      u[0] = 0.0 ;
   }
   else
   {
      fSecondDerivative[0] = -0.5 ;
      u[0] = (3.0/(fArgument[1]-fArgument[0]))*
	     ((fFunction[1]-fFunction[0])/(fArgument[1]-fArgument[0]) -
	                                                           pFirstDerStart) ;
   }
   
   // Decomposition loop for tridiagonal algorithm. fSecondDerivative[i] and u[i]
   // are used for temporary storage of the decomposed factors.
   
   for(i=1;i<fNumber-1;i++)
   {
      sig = (fArgument[i]-fArgument[i-1])/(fArgument[i+1]-fArgument[i-1]) ;
      p = sig*fSecondDerivative[i-1] + 2.0 ;
      fSecondDerivative[i] = (sig - 1.0)/p ;
      u[i] = (fFunction[i+1]-fFunction[i])/(fArgument[i+1]-fArgument[i]) -
	     (fFunction[i]-fFunction[i-1])/(fArgument[i]-fArgument[i-1]) ;
      u[i] =(6.0*u[i]/(fArgument[i+1]-fArgument[i-1]) - sig*u[i-1])/p ;
   }
   if(pFirstDerFinish > maxDerivative)
   {
      qn = 0.0 ;
      un = 0.0 ;
   }
   else
   {
      qn = 0.5 ;
      un =(3.0/(fArgument[fNumber-1]-fArgument[fNumber-2]))*(pFirstDerFinish -
	   (fFunction[fNumber-1]-fFunction[fNumber-2])/
	   (fArgument[fNumber-1]-fArgument[fNumber-2])) ;
   }
   fSecondDerivative[fNumber-1] = (un - qn*u[fNumber-2])/
                                  (qn*fSecondDerivative[fNumber-2] + 1.0) ;
   
   // The backsubstitution loop for the triagonal algorithm of solving a linear
   // system of equations.
   
   for(k=fNumber-2;k>=0;k--)
   {
      fSecondDerivative[k] = fSecondDerivative[k]*fSecondDerivative[k+1] + u[k] ;
   }
   delete[] u ;
} 

// ----------------------------------------------------------------------------
//
// Destructor deletes dynamically created arrays for data members: fArgument,
// fFunction and fSecondDerivative, all have dimension of fNumber
      
G4DataInterpolation::~G4DataInterpolation()
{
   delete[] fArgument ;
   delete[] fFunction ;
   if(fSecondDerivative) delete[] fSecondDerivative ;
}

// ------------------------------------------------------------------------
//
// This function returns the value P(pX), where P(x) is polynom of fNumber-1 degree
// such that P(fArgument[i]) = fFunction[i], for i = 0, ..., fNumber-1  . This is
// Lagrange's form of interpolation and it is based on Neville's algorithm

G4double 
G4DataInterpolation::PolynomInterpolation(G4double pX,
	                                  G4double& deltaY ) const
{
   G4int i, m, k = 0 ;
   G4double mult, diff, difi, deltaLow, deltaUp, cd, y ;
   G4double* c = new G4double[fNumber] ;
   G4double* d = new G4double[fNumber] ;
   diff = fabs(pX-fArgument[0]) ;
   for(i=1;i<fNumber;i++)
   {
      difi = fabs(pX-fArgument[i]) ;
      if(difi <diff)
      {
	 k = i ;
	 diff = difi ;
      }
      c[i] = fFunction[i] ;
      d[i] = fFunction[i] ;   
   }
   y = fFunction[k--] ;     
   for(m=1;m<fNumber;m++)
   {
      for(i=0;i<fNumber-m;i++)
      {
	 deltaLow = fArgument[i] - pX ;
	 deltaUp = fArgument[i+m] - pX ;
	 cd = c[i+1] - d[i] ;
	 mult = deltaLow - deltaUp ;
	 if(mult == 0.0)
	 {
	    G4Exception
	    ("Coincident nodes in G4DataInterpolation::PolynomInterpolation") ;
	 }
	 mult  = cd/mult ;
	 d[i] = deltaUp*mult ;
	 c[i] = deltaLow*mult ;
      }
      y += (deltaY = (2*k < (fNumber - m -1) ? c[k+1] : d[k--] )) ;
   }
   delete[] c ;
   delete[] d ;
   
   return y ;
}

// -----------------------------------------------------------
//
// Given arrays fArgument[0,..,fNumber-1] and fFunction[0,..,fNumber-1] , this
// function calculates an array of coefficients. The coefficients don't provide
// usually (fNumber>10) better accuracy for polynom interpolation, as compared with
// PolynomInterpolation function. They could be used instead for derivate 
// calculations and some other applications.

void 
G4DataInterpolation::PolIntCoefficient( G4double cof[]) const 
{
   G4int i, j ;
   G4double factor, reducedY, mult ;
   G4double* tempArgument = new G4double[fNumber] ;
   
   for(i=0;i<fNumber;i++)
   {
      tempArgument[i] = cof[i] = 0.0 ;
   }
   tempArgument[fNumber-1] = -fArgument[0] ;
   
   for(i=1;i<fNumber;i++)
   {
      for(j=fNumber-1-i;j<fNumber-1;j++)
      {
	 tempArgument[j] -= fArgument[i]*tempArgument[j+1] ;
      }
      tempArgument[fNumber-1] -= fArgument[i] ;
   }
   for(i=0;i<fNumber;i++)
   {
      factor = fNumber ;
      for(j=fNumber-1;j>=1;j--)
      {
	 factor = j*tempArgument[j] + factor*fArgument[i] ;
      }
      reducedY = fFunction[i]/factor ;
      mult = 1.0 ;
      for(j=fNumber-1;j>=0;j--)
      {
	 cof[j] += mult*reducedY ;
	 mult = tempArgument[j] + mult*fArgument[i] ;
      }
   }
   delete[] tempArgument ;
}



// ----------------------------------------------------------------

// The function returns diagonal rational function (Bulirsch and Stoer algorithm
// of Neville type) Pn(x)/Qm(x) where P and Q are polynoms.
// Tests showed the method is not stable and hasn't advantage if compared with
// polynomial interpolation ?!


G4double
G4DataInterpolation::RationalPolInterpolation(G4double pX,
	                                      G4double& deltaY ) const 
{
   G4int i, m, k = 0 ;
   const G4double tolerance = 1.6e-24 ;
   G4double mult, difi, diff, cd, y, cof ;
   G4double* c = new G4double[fNumber] ;
   G4double* d = new G4double[fNumber] ;
   diff = fabs(pX-fArgument[0]) ;
   for(i=0;i<fNumber;i++)
   {
      difi = fabs(pX-fArgument[i]) ;
      if(difi == 0.0)
      {
	 y = fFunction[i] ;
	 deltaY = 0.0 ;
         delete[] c ;
         delete[] d ;
	 return y ;
      }
      else if(difi < diff)
      {
	 k = i ;
	 diff = difi ;
      }
      c[i] = fFunction[i] ;
      d[i] = fFunction[i] + tolerance ;   // to prevent rare zero/zero cases
   }
   y = fFunction[k--] ;    
   for(m=1;m<fNumber;m++)
   {
      for(i=0;i<fNumber-m;i++)
      {
	 cd  = c[i+1] - d[i] ;
	 difi  = fArgument[i+m] - pX ;
	 cof  = (fArgument[i] - pX)*d[i]/difi ;
	 mult = cof - c[i+1] ;
	 if(mult == 0.0)    // function to be interpolated has pole at pX
	 {
	    G4Exception("Error in G4DataInterpolation::RationalPolInterpolation") ;
	 }
	 mult  = cd/mult ;
	 d[i] = c[i+1]*mult ;
	 c[i] = cof*mult ;
      }
      y += (deltaY = (2*k < (fNumber - m - 1) ? c[k+1] : d[k--] )) ;
   }
   delete[] c ;
   delete[] d ;
   
   return y ;
}

// --------------------------------------------------------------------------

// Cubic spline interpolation in point pX for function given by the table:
// fArgument, fFunction. The constructor, which creates fSecondDerivative, must be
// called before. The function works optimal, if sequential calls are in random
// values of pX.

G4double 
G4DataInterpolation::CubicSplineInterpolation(G4double pX) const 
{
   G4int kLow, kHigh, k ;
   G4double deltaHL, a, b ;
   
   // Searching in the table by means of bisection method. 
   // fArgument must be monotonic, either increasing or decreasing
   
   kLow = 0 ;
   kHigh = fNumber - 1 ;
   while((kHigh - kLow) > 1)
   {
      k = (kHigh + kLow) >> 1 ; // compute midpoint 'bisection'
      if(fArgument[k] > pX)
      {
	 kHigh = k ;
      }
      else
      {
	 kLow = k ;
      }
   }        // kLow and kHigh now bracket the input value of pX
   deltaHL = fArgument[kHigh] - fArgument[kLow] ;
   if(deltaHL == 0.0)
   {
      G4Exception(
      "Bad fArgument input in G4DataInterpolation::CubicSplineInterpolation") ;
   }
   a = (fArgument[kHigh] - pX)/deltaHL ;
   b = (pX - fArgument[kLow])/deltaHL ;
   
   // Final evaluation of cubic spline polynomial for return
   
   return a*fFunction[kLow] + b*fFunction[kHigh] + 
          ((a*a*a - a)*fSecondDerivative[kLow] +
           (b*b*b - b)*fSecondDerivative[kHigh])*deltaHL*deltaHL/6.0  ;
}

// ---------------------------------------------------------------------
//
// Return cubic spline interpolation in the point pX which is located between
// fArgument[index] and fArgument[index+1]. It is usually called in sequence of
// known from external analysis values of index.


G4double 
G4DataInterpolation::FastCubicSpline(G4double pX, 
				     G4int index) const 
{
   G4double delta, a, b ;
   delta = fArgument[index+1] - fArgument[index] ;
   if(delta == 0.0)
   {
      G4Exception("Bad fArgument input in G4DataInterpolation::FastCubicSpline") ;
   }
   a = (fArgument[index+1] - pX)/delta ;
   b = (pX - fArgument[index])/delta ;
   
   // Final evaluation of cubic spline polynomial for return
   
   return a*fFunction[index] + b*fFunction[index+1] + 
          ((a*a*a - a)*fSecondDerivative[index] +
           (b*b*b - b)*fSecondDerivative[index+1])*delta*delta/6.0  ;
}

// ---------------------------------------------------------------------------
//
// Given argument pX, returns index k, so that pX bracketed by fArgument[k] and
// fArgument[k+1]

G4int 
G4DataInterpolation::LocateArgument(G4double pX) const 
{
   G4int kLow, kHigh, k ; 
   G4bool ascend ;
   kLow = -1 ;
   kHigh = fNumber ;
   ascend = (fArgument[fNumber-1] >= fArgument[0]) ;
   while((kHigh - kLow) > 1)
   {
      k = (kHigh + kLow) >> 1 ; // compute midpoint 'bisection'
      if(pX >= fArgument[k] == ascend)
      {
	 kLow = k ;
      }
      else
      {
	 kHigh = k ;
      }
   }
   if(pX == fArgument[0])
   {
      return 1 ;
   }
   else if(pX == fArgument[fNumber-1]) 
   {
      return fNumber - 2 ;
   }
   else  return kLow ;
}

// ------------------------------------------------------------------------
//
// Given a value pX, returns a value 'index' such that pX is between fArgument[index]
// and fArgument[index+1]. fArgument MUST BE MONOTONIC, either increasing or
// decreasing. If index = -1 or fNumber, this indicates that pX is out of range.
// The value index on input is taken as the initial approximation for index on
// output.


void 
G4DataInterpolation::CorrelatedSearch( G4double pX,
	                               G4int& index ) const 
{
   G4int kHigh, k, Increment ;
   // ascend = true for ascending order of table, false otherwise
   G4bool ascend = (fArgument[fNumber-1] >= fArgument[0]) ;
   if(index < 0 || index > fNumber-1)
   {
      index = -1 ;
      kHigh = fNumber ;
   }
   else
   {
      Increment = 1 ;   //   What value would be the best ?
      if((pX >= fArgument[index]) == ascend)
      {
	 if(index == fNumber -1)
	 {
	    index = fNumber ;
	    return ;
	 }
	 kHigh = index + 1 ;
	 while((pX >= fArgument[kHigh]) == ascend)
	 {
	    index = kHigh ;
	    Increment += Increment ;      // double the Increment
	    kHigh = index + Increment ;
	    if(kHigh > (fNumber - 1))
	    {
	       kHigh = fNumber ;
	       break ;
	    }
	 }
      }
      else
      {
	 if(index == 0)
	 {
	    index = -1 ;
	    return ;
	 }
         kHigh = index-- ;
         while((pX < fArgument[index]) == ascend)
         {
	    kHigh = index ;
	    Increment <<= 1 ;      // double the Increment
	    if(Increment >= kHigh)
	    {
	       index = -1 ;
	       break ;
	    }
	    else
	    {
	       index = kHigh - Increment ;
	    }
         }
      }        // Value bracketed
   }
   // final bisection searching

   while((kHigh - index) != 1)
   {
      k = (kHigh + index) >> 1 ;
      if((pX >= fArgument[k]) == ascend)
      {
	 index = k ;
      }
      else
      {
	 kHigh = k ;
      }
   }
   if(pX == fArgument[fNumber-1])
   {
      index = fNumber - 2 ;
   }
   if(pX == fArgument[0])
   {
      index = 0 ;
   }
   return ;
}
