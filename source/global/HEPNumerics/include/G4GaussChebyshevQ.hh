// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GaussChebyshevQ.hh,v 1.3 2000-11-20 17:26:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// Class for Gauss-Chebyshev quadrature method
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// ------------------------------ CONSTRUCTORS ----------------------------
//
// Constructor for Gauss-Chebyshev quadrature method
//
// G4GaussChebyshevQuadrature( function pFunction,
//			       G4int nChebyshev       )
//
//
//
// ------------------------------- METHODS -----------------------------------
//
// Integrates function pointed by fFunction from a to b by Gauss-Chebyshev quadrature
// method
//
// G4double Integral(G4double a, G4double b) const 

// ------------------------------- HISTORY --------------------------------
//
// 13.05.97 V.Grichine (Vladimir.Grichine@cern.ch)

#ifndef G4GAUSSCHEBYSHEVQ_HH
#define G4GAUSSCHEBYSHEVQ_HH

#include "G4VGaussianQuadrature.hh"

class G4GaussChebyshevQ : public G4VGaussianQuadrature
{
public:
        // Constructor/destructor

        G4GaussChebyshevQ( function pFunction,
			   G4int nChebyshev           ) ;
			   
	~G4GaussChebyshevQ() ;		   
			       
        // Methods
			     
        G4double Integral(G4double a, G4double b) const ;


private:

	G4GaussChebyshevQ(const G4GaussChebyshevQ&);
	G4GaussChebyshevQ& operator=(const G4GaussChebyshevQ&);

};

#endif
