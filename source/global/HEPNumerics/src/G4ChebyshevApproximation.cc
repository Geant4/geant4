//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4ChebyshevApproximation.cc 67970 2013-03-13 10:10:06Z gcosmo $
//

#include "G4ChebyshevApproximation.hh"
#include "G4PhysicalConstants.hh"

// Constructor for initialisation of the class data members.
// It creates the array fChebyshevCof[0,...,fNumber-1], fNumber = n ;
// which consists of Chebyshev coefficients describing the function
// pointed by pFunction. The values a and b fix the interval of validity
// of the Chebyshev approximation. 

G4ChebyshevApproximation::G4ChebyshevApproximation( function pFunction,
                                                    G4int n, 
                                                    G4double a,
                                                    G4double b       ) 
   : fFunction(pFunction), fNumber(n),
     fChebyshevCof(new G4double[fNumber]),
     fMean(0.5*(b+a)), fDiff(0.5*(b-a))
{
   G4int i=0, j=0 ;
   G4double  rootSum=0.0, cofj=0.0 ;
   G4double* tempFunction = new G4double[fNumber] ;
   G4double weight = 2.0/fNumber ;
   G4double cof = 0.5*weight*pi ;    // pi/n
   
   for (i=0;i<fNumber;i++)
   {
      rootSum = std::cos(cof*(i+0.5)) ;
      tempFunction[i]= fFunction(rootSum*fDiff+fMean) ;
   }
   for (j=0;j<fNumber;j++) 
   {
      cofj = cof*j ;
      rootSum = 0.0 ;
      
      for (i=0;i<fNumber;i++)
      {
         rootSum += tempFunction[i]*std::cos(cofj*(i+0.5)) ;
      }
      fChebyshevCof[j] = weight*rootSum ;
   }
   delete[] tempFunction ;
}

// --------------------------------------------------------------------
//
// Constructor for creation of Chebyshev coefficients for mx-derivative
// from pFunction. The value of mx ! MUST BE ! < nx , because the result
// array of fChebyshevCof will be of (nx-mx) size.  The values a and b
// fix the interval of validity of the Chebyshev approximation. 

G4ChebyshevApproximation::
G4ChebyshevApproximation( function pFunction,
                          G4int nx, G4int mx,
                          G4double a, G4double b ) 
   : fFunction(pFunction), fNumber(nx),
     fChebyshevCof(new G4double[fNumber]),
     fMean(0.5*(b+a)), fDiff(0.5*(b-a))
{
   if(nx <= mx)
   {
      G4Exception("G4ChebyshevApproximation::G4ChebyshevApproximation()",
                  "InvalidCall", FatalException, "Invalid arguments !") ;
   }
   G4int i=0, j=0 ;
   G4double  rootSum = 0.0, cofj=0.0;   
   G4double* tempFunction = new G4double[fNumber] ;
   G4double weight = 2.0/fNumber ;
   G4double cof = 0.5*weight*pi ;    // pi/nx
   
   for (i=0;i<fNumber;i++)
   {
      rootSum = std::cos(cof*(i+0.5)) ;
      tempFunction[i] = fFunction(rootSum*fDiff+fMean) ;
   }
   for (j=0;j<fNumber;j++) 
   {
      cofj = cof*j ;
      rootSum = 0.0 ;
      
      for (i=0;i<fNumber;i++)
      {
         rootSum += tempFunction[i]*std::cos(cofj*(i+0.5)) ;
      }
      fChebyshevCof[j] = weight*rootSum ; // corresponds to pFunction
   }
   // Chebyshev coefficients for (mx)-derivative of pFunction
   
   for(i=1;i<=mx;i++)
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
   : fFunction(pFunction), fNumber(n),
     fChebyshevCof(new G4double[fNumber]),
     fMean(0.5*(b+a)), fDiff(0.5*(b-a))
{
   G4int i=0, j=0;
   G4double  rootSum=0.0, cofj=0.0;
   G4double* tempFunction = new G4double[fNumber] ;
   G4double weight = 2.0/fNumber;
   G4double cof = 0.5*weight*pi ;    // pi/n
   
   for (i=0;i<fNumber;i++)
   {
      rootSum = std::cos(cof*(i+0.5)) ;
      tempFunction[i]= fFunction(rootSum*fDiff+fMean) ;
   }
   for (j=0;j<fNumber;j++) 
   {
      cofj = cof*j ;
      rootSum = 0.0 ;
      
      for (i=0;i<fNumber;i++)
      {
         rootSum += tempFunction[i]*std::cos(cofj*(i+0.5)) ;
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
      G4Exception("G4ChebyshevApproximation::GetChebyshevCof()",
                  "InvalidCall", FatalException, "Argument out of range !") ;
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
   G4double evaluate = 0.0, evaluate2 = 0.0, temp = 0.0,
            xReduced = 0.0, xReduced2 = 0.0 ;

   if ((x-fMean+fDiff)*(x-fMean-fDiff) > 0.0) 
   {
      G4Exception("G4ChebyshevApproximation::ChebyshevEvaluation()",
                  "InvalidCall", FatalException, "Invalid argument !") ;
   }
   xReduced = (x-fMean)/fDiff ;
   xReduced2 = 2.0*xReduced ;
   for (G4int i=fNumber-1;i>=1;i--) 
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
   G4double cof = 1.0/fDiff ;
   derCof[fNumber-1] = 0.0 ;
   derCof[fNumber-2] = 2*(fNumber-1)*fChebyshevCof[fNumber-1] ;
   for(G4int i=fNumber-3;i>=0;i--)
   {
      derCof[i] = derCof[i+2] + 2*(i+1)*fChebyshevCof[i+1] ;  
   }
   for(G4int j=0;j<fNumber;j++)
   {
      derCof[j] *= cof ;
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
   G4double cof = 0.5*fDiff, sum = 0.0, factor = 1.0 ;
   for(G4int i=1;i<fNumber-1;i++)
   {
      integralCof[i] = cof*(fChebyshevCof[i-1] - fChebyshevCof[i+1])/i ;
      sum += factor*integralCof[i] ;
      factor = -factor ;
   }
   integralCof[fNumber-1] = cof*fChebyshevCof[fNumber-2]/(fNumber-1) ;
   sum += factor*integralCof[fNumber-1] ;
   integralCof[0] = 2.0*sum ;                // set the constant of integration
}                                         
