// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ChebyshevApproximation.cc,v 1.1 1999-01-07 16:08:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ChebyshevApproximation.hh"


// Constructor for initialisation of the class data members. It creates the array
// fChebyshevCof[0,...,fNumber-1], fNumber = n ; which consists of Chebyshev
// coefficients describing the function pointed by pFunction. The values a and b
// fixe the interval of validity of Chebyshev approximation. 


G4ChebyshevApproximation::G4ChebyshevApproximation( function pFunction,
                                                    G4int n, 
                                                    G4double a,
			                            G4double b       ) 
{
   G4int i, j ;
   G4double  rootSum, cof, cofj, weight ;
   
   fFunction = pFunction ;
   fNumber = n ;
   fDiff = 0.5*(b-a) ;
   fMean = 0.5*(b+a) ;
   fChebyshevCof = new G4double[fNumber] ;
   
   G4double* tempFunction = new G4double[fNumber] ;
   
   weight = 2.0/fNumber ;
   cof = 0.5*weight*pi ;    // pi/n
   
   for (i=0;i<fNumber;i++)
   {
      rootSum = cos(cof*(i+0.5)) ;
      tempFunction[i]= fFunction(rootSum*fDiff+fMean) ;
   }
   for (j=0;j<fNumber;j++) 
   {
      cofj = cof*j ;
      rootSum = 0.0 ;
      
      for (i=0;i<fNumber;i++)
      {
         rootSum += tempFunction[i]*cos(cofj*(i+0.5)) ;
      }
      fChebyshevCof[j] = weight*rootSum ;
   }
   delete[] tempFunction ;
}

// --------------------------------------------------------------------
//
// Constructor for creation of Chebyshev coefficients for m-derivative
// from pFunction. The value of m ! MUST BE ! < n , because the result
// array of fChebyshevCof will be of (n-m) size.  The values a and b
// fixe the interval of validity of Chebyshev approximation. 

	
G4ChebyshevApproximation::
G4ChebyshevApproximation( function pFunction,
                          G4int n,
			  G4int m,
                          G4double a,
			  G4double b       ) 
{
   if(n <= m)
   {
      G4Exception
      ("Invalid arguments in G4ChebyshevApproximation::G4ChebyshevApproximation") ;
   }
   G4int i, j ;
   G4double  rootSum, cof, cofj, weight ;
   
   fFunction = pFunction ;
   fNumber = n ;
   fDiff = 0.5*(b-a) ;
   fMean = 0.5*(b+a) ;
   fChebyshevCof = new G4double[fNumber] ;
   
   G4double* tempFunction = new G4double[fNumber] ;
   
   weight = 2.0/fNumber ;
   cof = 0.5*weight*pi ;    // pi/n
   
   for (i=0;i<fNumber;i++)
   {
      rootSum = cos(cof*(i+0.5)) ;
      tempFunction[i] = fFunction(rootSum*fDiff+fMean) ;
   }
   for (j=0;j<fNumber;j++) 
   {
      cofj = cof*j ;
      rootSum = 0.0 ;
      
      for (i=0;i<fNumber;i++)
      {
         rootSum += tempFunction[i]*cos(cofj*(i+0.5)) ;
      }
      fChebyshevCof[j] = weight*rootSum ; // corresponds to pFunction
   }
   // Chebyshev coefficients for (m)-derivative of pFunction
   
   for(i=1;i<=m;i++)
   {
      DerivativeChebyshevCof(tempFunction) ;
      fNumber-- ;
      for(j=0;j<fNumber;j++)
      {
	 fChebyshevCof[j] = tempFunction[j] ; // corresponds to (i)-derivative
      }
   }
   delete[] tempFunction ;   // delete of dynamically allocated tempFunction
}

// ------------------------------------------------------
//
// Constructor for creation of Chebyshev coefficients for integral
// from pFunction.
	
G4ChebyshevApproximation::G4ChebyshevApproximation( function pFunction,
                                                    G4double a,
			                            G4double b, 
                                                    G4int n            ) 
{
   G4int i,j;
   G4double  rootSum, cof, cofj, weight ;
   
   fFunction = pFunction ;
   fNumber = n ;
   fDiff = 0.5*(b-a) ;
   fMean = 0.5*(b+a) ;
   fChebyshevCof = new G4double[fNumber] ;
   
   G4double* tempFunction = new G4double[fNumber] ;
   
   weight = 2.0/fNumber ;
   cof = 0.5*weight*pi ;    // pi/n
   
   for (i=0;i<fNumber;i++)
   {
      rootSum = cos(cof*(i+0.5)) ;
      tempFunction[i]= fFunction(rootSum*fDiff+fMean) ;
   }
   for (j=0;j<fNumber;j++) 
   {
      cofj = cof*j ;
      rootSum = 0.0 ;
      
      for (i=0;i<fNumber;i++)
      {
         rootSum += tempFunction[i]*cos(cofj*(i+0.5)) ;
      }
      fChebyshevCof[j] = weight*rootSum ; // corresponds to pFunction
   }
   // Chebyshev coefficients for integral of pFunction
   
   IntegralChebyshevCof(tempFunction) ;
   for(j=0;j<fNumber;j++)
   {
      fChebyshevCof[j] = tempFunction[j] ; // corresponds to integral
   }
   delete[] tempFunction ;   // delete of dynamically allocated tempFunction
}



// ---------------------------------------------------------------
//
// Destructor deletes the array of Chebyshev coefficients

G4ChebyshevApproximation::~G4ChebyshevApproximation()
{
   delete[] fChebyshevCof ;
}

// ---------------------------------------------------------------
//
// Access function for Chebyshev coefficients
//


G4double
G4ChebyshevApproximation::GetChebyshevCof(G4int number) const 
{
   if(number < 0 && number >= fNumber)
   {
      G4Exception
      ("Argument out of range in G4ChebyshevApproximation::GetChebyshevCof") ;
   }
   return fChebyshevCof[number] ;
}

// --------------------------------------------------------------
//
// Evaluate the value of fFunction at the point x via the Chebyshev coefficients
// fChebyshevCof[0,...,fNumber-1]

G4double
G4ChebyshevApproximation::ChebyshevEvaluation(G4double x) const 
{
	G4int i;
	G4double evaluate = 0.0, evaluate2 = 0.0, temp, xReduced, xReduced2 ;

	if ((x-fMean+fDiff)*(x-fMean-fDiff) > 0.0) 
	{
G4Exception("Invalid argument in G4ChebyshevApproximation::ChebyshevEvaluation");
	}
	xReduced = (x-fMean)/fDiff ;
	xReduced2 = 2.0*xReduced ;
	for (i=fNumber-1;i>=1;i--) 
	{
	   temp = evaluate ;
	   evaluate  = xReduced2*evaluate - evaluate2 + fChebyshevCof[i] ;
	   evaluate2 = temp ;
	}
	return xReduced*evaluate - evaluate2 + 0.5*fChebyshevCof[0] ;
}

// ------------------------------------------------------------------
//
// Returns the array derCof[0,...,fNumber-2], the Chebyshev coefficients of the 
// derivative of the function whose coefficients are fChebyshevCof

void
G4ChebyshevApproximation::DerivativeChebyshevCof(G4double derCof[]) const 
{
   G4int i ;
   G4double cof = 1.0/fDiff ;
   derCof[fNumber-1] = 0.0 ;
   derCof[fNumber-2] = 2*(fNumber-1)*fChebyshevCof[fNumber-1] ;
   for(i=fNumber-3;i>=0;i--)
   {
      derCof[i] = derCof[i+2] + 2*(i+1)*fChebyshevCof[i+1] ;  
   }
   for(i=0;i<fNumber;i++)
   {
      derCof[i] *= cof ;
   }
}

// ------------------------------------------------------------------------
//
// This function produces the array integralCof[0,...,fNumber-1] , the Chebyshev
// coefficients of the integral of the function whose coefficients are  
// fChebyshevCof[]. The constant of integration is set so that the integral 
// vanishes at the point (fMean - fDiff), i.e. at the begining of the interval of
// validity (we start the integration from this point).
//
   
void 
G4ChebyshevApproximation::IntegralChebyshevCof(G4double integralCof[]) const 
{
   G4int i ;
   G4double cof = 0.5*fDiff, sum = 0.0, factor = 1.0 ;
   for(i=1;i<fNumber-1;i++)
   {
      integralCof[i] = cof*(fChebyshevCof[i-1] - fChebyshevCof[i+1])/i ;
      sum += factor*integralCof[i] ;
      factor = -factor ;
   }
   integralCof[fNumber-1] = cof*fChebyshevCof[fNumber-2]/(fNumber-1) ;
   sum += factor*integralCof[fNumber-1] ;
   integralCof[0] = 2.0*sum ;                // set the constant of integration
}                                         
