// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GaussLaguerreQ.cc,v 1.3 2000-11-20 17:26:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4GaussLaguerreQ.hh"



// ------------------------------------------------------------
//
// Constructor for Gauss-Laguerre quadrature method: integral from zero to
// infinity of pow(x,alpha)*exp(-x)*f(x). The value of nLaguerre sets the accuracy.
// The constructor creates arrays fAbscissa[0,..,nLaguerre-1] and 
// fWeight[0,..,nLaguerre-1] . 
//

G4GaussLaguerreQ::G4GaussLaguerreQ( function pFunction,
				    G4double alpha,
			            G4int nLaguerre      ) 
   : G4VGaussianQuadrature(pFunction)
{
   const G4double tolerance = 1.0e-10 ;
   const G4int maxNumber = 12 ;
   G4int i, j, k ;
   G4double newton=0.;
   G4double newton1, temp1, temp2, temp3, temp, cofi ;

   fNumber = nLaguerre ;
   fAbscissa = new G4double[fNumber] ;
   fWeight   = new G4double[fNumber] ;
      
   for(i=1;i<=fNumber;i++)      // Loop over the desired roots
   {
      if(i == 1)
      {
	 newton = (1.0 + alpha)*(3.0 + 0.92*alpha)/(1.0 + 2.4*fNumber + 1.8*alpha) ;
      }
      else if(i == 2)
      {
	 newton += (15.0 + 6.25*alpha)/(1.0 + 0.9*alpha + 2.5*fNumber) ;
      }
      else
      {
	 cofi = i - 2 ;
	 newton += ((1.0+2.55*cofi)/(1.9*cofi) + 1.26*cofi*alpha/(1.0+3.5*cofi))*
	           (newton - fAbscissa[i-3])/(1.0 + 0.3*alpha) ;
      }
      for(k=1;k<=maxNumber;k++)
      {
	 temp1 = 1.0 ;
	 temp2 = 0.0 ;
	 for(j=1;j<=fNumber;j++)
	 {
	    temp3 = temp2 ;
	    temp2 = temp1 ;
	    temp1 = ((2*j - 1 + alpha - newton)*temp2 - (j - 1 + alpha)*temp3)/j ;
	 }
	 temp = (fNumber*temp1 - (fNumber +alpha)*temp2)/newton ;
	 newton1 = newton ;
	 newton  = newton1 - temp1/temp ;
         if(fabs(newton - newton1) <= tolerance) 
	 {
	    break ;
	 }
      }
      if(k > maxNumber)
      {
	 G4Exception("Too many iterations in Gauss-Laguerre constructor") ;
      }
	 
      fAbscissa[i-1] =  newton ;
      fWeight[i-1] = -exp(GammaLogarithm(alpha + fNumber) - 
			  GammaLogarithm((G4double)fNumber))/(temp*fNumber*temp2) ;
   }
}

// -----------------------------------------------------------------
//
// Gauss-Laguerre method for integration of pow(x,alpha)*exp(-x)*pFunction(x)
// from zero up to infinity. pFunction is evaluated in fNumber points for which
// fAbscissa[i] and fWeight[i] arrays were created in
// G4VGaussianQuadrature(double,int) constructor

G4double 
G4GaussLaguerreQ::Integral() const 
{
   G4int i ;
   G4double integral = 0.0 ;
   for(i=0;i<fNumber;i++)
   {
      integral += fWeight[i]*fFunction(fAbscissa[i]) ;
   }
   return integral ;
}
