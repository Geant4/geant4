// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ChebyshevApproximation.hh,v 1.1 1999-01-07 16:08:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class creating the Chebyshev approximation for a function pointed by fFunction
// data member. The Chebyshev polinom approximation provides an efficient evaluation
// of minimax polynomial, which (among all polynomials of the same degree) has the
// smallest maximum deviation from the true function. 
// The methods based mainly on recommendations given in the book : An introduction to
// NUMERICAL METHODS IN C++, B.H. Flowers, Claredon Press, Oxford, 1995
//
// ------------------------- MEMBER DATA ------------------------------------
//
// function   fFunction - pointer to a function considered
// G4int      fNumber - number of Chebyshev coefficients
// G4double*  fChebyshevCof - array of Chebyshev coefficients
// G4double   fMean = (a+b)/2 - mean point of interval
// G4double   fDiff = (b-a)/2 - half of the interval value
//
// ------------------------ CONSTRUCTORS ----------------------------------
//
// Constructor for initialisation of the class data members. It creates the array
// fChebyshevCof[0,...,fNumber-1], fNumber = n ; which consists of Chebyshev
// coefficients describing the function pointed by pFunction. The values a and b
// fixe the interval of validity of Chebyshev approximation.
//
// G4ChebyshevApproximation( function pFunction,
//                           G4int n, 
//                           G4double a,
// 			     G4double b       ) 
// 	
// --------------------------------------------------------------------
//
// Constructor for creation of Chebyshev coefficients for m-derivative
// from pFunction. The value of m ! MUST BE ! < n , because the result
// array of fChebyshevCof will be of (n-m) size. There is a definite dependence
// between the proper selection of n, m, a and b values to get better accuracy
// of the derivative value.
//	
// G4ChebyshevApproximation( function pFunction,
//                           G4int n,
// 			     G4int m,
//                           G4double a,
//			     G4double b       ) 
//
// ------------------------------------------------------
//
// Constructor for creation of Chebyshev coefficients for integral
// from pFunction.
//	
// G4ChebyshevApproximation( function pFunction,
//                           G4double a,
//			     G4double b, 
//                           G4int n            ) 
//
// ---------------------------------------------------------------
//
// Destructor deletes the array of Chebyshev coefficients
//
// ~G4ChebyshevApproximation()
//
// ----------------------------- METHODS ----------------------------------
//
// Access function for Chebyshev coefficients
//
// G4double GetChebyshevCof(G4int number) const 
//
// --------------------------------------------------------------
//
// Evaluate the value of fFunction at the point x via the Chebyshev coefficients
// fChebyshevCof[0,...,fNumber-1]
//
// G4double ChebyshevEvaluation(G4double x) const 
//
// ------------------------------------------------------------------
//
// Returns the array derCof[0,...,fNumber-2], the Chebyshev coefficients of the 
// derivative of the function whose coefficients are fChebyshevCof
//
// void DerivativeChebyshevCof(G4double derCof[]) const 
//
// ------------------------------------------------------------------------
//
// This function produces the array integralCof[0,...,fNumber-1] , the Chebyshev
// coefficients of the integral of the function whose coefficients are 
// fChebyshevCof. The constant of integration is set so that the integral vanishes
// at the point (fMean - fDiff)
//   
// void IntegralChebyshevCof(G4double integralCof[]) const 
//
// --------------------------- HISTORY --------------------------------------
//
//  24.04.97   V.Grichine ( Vladimir.Grichine@cern.ch )

#ifndef G4CHEBYSHEVAPPROXIMATION_HH
#define G4CHEBYSHEVAPPROXIMATION_HH

#include "globals.hh"

typedef G4double (*function)(G4double) ;

class G4ChebyshevApproximation
{
public:
        G4ChebyshevApproximation( function pFunction,
                                  G4int n, 
                                  G4double a,
			          G4double b       ) ;
        
	// Constructor for creation of Chebyshev coefficients for m-derivative
	// from pFunction. The value of m ! MUST BE ! < n , because the result
	// array of fChebyshevCof will be of (n-m) size.
	
	G4ChebyshevApproximation( function pFunction,
                                  G4int n,
				  G4int m,
                                  G4double a,
			          G4double b       ) ;

	// Constructor for creation of Chebyshev coefficients for integral
	// from pFunction.
	
	G4ChebyshevApproximation( function pFunction,
                                  G4double a,
			          G4double b, 
                                  G4int n            ) ;
				  
				  
       
       ~G4ChebyshevApproximation() ;
       
        // Access functions
       
        G4double GetChebyshevCof(G4int number) const ;
       
        // Methods
		
	G4double ChebyshevEvaluation(G4double x) const ;	
	
	void DerivativeChebyshevCof(G4double derCof[]) const ;
	
	void IntegralChebyshevCof(G4double integralCof[]) const ;
       
protected:

private:

        function   fFunction ;
	G4int      fNumber ;
	G4double*  fChebyshevCof ;
	G4double   fMean ;
	G4double   fDiff ;

} ;

#endif
