// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGaussianQuadrature.hh,v 1.3 2000-11-20 17:26:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// Base Class for realisation of numerical methodes for integration of functions
// with signature double f(double) by Gaussian quadrature methods
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// ---------------------------- Member data: ----------------------------------
//
//  fFunction  - pointer to the function to be integrated
//  fNumber    - the number of points in fAbscissa and fWeight arrays
//  fAbscissa  - array of abscissas, where function will be evaluated
//  fWeight    - array of corresponding weights
//
//
// ----------------------------------------------------------------------
//
// Auxiliary function which returns the value of log(gamma-function(x))
//
// G4double 
// GammaLogarithm(G4double xx)

// ------------------------------------------------------------------------------
//
// History:
//             18.04.97   V.Grichine ( Vladimir.Grichine@cern.ch )

#ifndef G4VGAUSSIANQUADRATURE_HH
#define G4VGAUSSIANQUADRATURE_HH

#include "globals.hh"

typedef G4double (*function)(G4double) ;

class G4VGaussianQuadrature
{
public:
           // Base constructor

           G4VGaussianQuadrature( function pFunction ) ;
       
	   // Virtual destructor		     
			     
           virtual ~G4VGaussianQuadrature() ;
      
           // Access functions:
       
           G4double GetAbscissa(G4int index) const ;

           G4double GetWeight(G4int index) const ;
	   
	   G4int GetNumber() const { return fNumber ; }
       
           // Methods:
       
           // virtual G4double DefiniteIntegral( G4double a,
	   //                                    G4double b   ) const = 0 ;
       
           // virtual G4double Integral() const = 0 ;
       
 
 			    
protected:
           G4double GammaLogarithm(G4double xx) ;

           //  Data members common for GaussianQuadrature family
	   
	   function  fFunction ;
	   G4double* fAbscissa ;
	   G4double* fWeight ;
	   G4int     fNumber ;
private:

	   G4VGaussianQuadrature(const G4VGaussianQuadrature&);
	   G4VGaussianQuadrature& operator=(const G4VGaussianQuadrature&);

};

#endif
