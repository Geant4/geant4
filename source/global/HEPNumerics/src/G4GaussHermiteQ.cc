// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GaussHermiteQ.cc,v 1.2 1999-11-16 17:31:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4GaussHermiteQ.hh"


// ----------------------------------------------------------
//
// Constructor for Gauss-Hermite

G4GaussHermiteQ::G4GaussHermiteQ(      function pFunction, 
				       G4int nHermite       ) 
   : G4VGaussianQuadrature(pFunction)
{
   const G4double tolerance = 1.0e-12 ;
   const G4int maxNumber = 12 ;
   
   G4int i, j, k ;
   G4double newton, newton1, temp1, temp2, temp3, temp ;
   G4double piInMinusQ = pow(pi,-0.25) ;                // 1.0/sqrt(sqrt(pi)) ??

   fNumber = (nHermite +1)/2 ;
   fAbscissa = new G4double[fNumber] ;
   fWeight   = new G4double[fNumber] ;

   for(i=1;i<=fNumber;i++)
   {
      if(i == 1)
      {
	 newton = sqrt((G4double)(2*nHermite + 1)) - 
	          1.85575001*pow((G4double)(2*nHermite + 1),-0.16666999) ;
      }
      else if(i == 2)
      {
	 newton -= 1.14001*pow((G4double)nHermite,0.425999)/newton ;
      }
      else if(i == 3)
      {
	 newton = 1.86002*newton - 0.86002*fAbscissa[0] ;
      }
      else if(i == 4)
      {
	 newton = 1.91001*newton - 0.91001*fAbscissa[1] ;
      }
      else 
      {
	 newton = 2.0*newton - fAbscissa[i - 3] ;
      }
      for(k=1;k<=maxNumber;k++)
      {
	 temp1 = piInMinusQ ;
	 temp2 = 0.0 ;
	 for(j=1;j<=nHermite;j++)
	 {
	    temp3 = temp2 ;
	    temp2 = temp1 ;
	    temp1 = newton*sqrt(2.0/j)*temp2 - sqrt(((G4double)(j - 1))/j)*temp3 ;
	 }
	 temp = sqrt((G4double)2*nHermite)*temp2 ;
	 newton1 = newton ;
	 newton = newton1 - temp1/temp ;
         if(fabs(newton - newton1) <= tolerance) 
	 {
	    break ;
	 }
      }
      if(k > maxNumber)
      {
	 G4Exception("Too many iterations in Gauss-Hermite constructor") ;
      }
      fAbscissa[i-1] =  newton ;
      fWeight[i-1] = 2.0/(temp*temp) ;
   }
}


// ----------------------------------------------------------
//
// Gauss-Hermite method for integration of exp(-x*x)*nFunction(x) from minus infinity
// to plus infinity . 

G4double 
   G4GaussHermiteQ::Integral() const 
{
   G4int i ;
   G4double integral = 0.0 ;
   for(i=0;i<fNumber;i++)
   {
      integral += fWeight[i]*(fFunction(fAbscissa[i]) + fFunction(-fAbscissa[i])) ;
   }
   return integral ;
}
