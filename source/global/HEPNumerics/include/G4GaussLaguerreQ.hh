// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GaussLaguerreQ.hh,v 1.2 1999-11-16 17:30:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// Class for realization of Gauss-Laguerre quadrature method
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book: M. Abramowitz, I. Stegun, Handbook
// of mathematical functions, DOVER Publications INC, New York 1965 ; chapters 9,
// 10, and 22 .
//
// ---------------------------------------------------------------------------
//
// Constructor for Gauss-Laguerre quadrature method: integral from zero to
// infinity of pow(x,alpha)*exp(-x)*f(x). The value of nLaguerre sets the accuracy.
// The constructor creates arrays fAbscissa[0,..,nLaguerre-1] and 
// fWeight[0,..,nLaguerre-1] . The function GaussLaguerre(f) should be called
// then with any f .
//
// G4GaussLaguerreQ( function pFunction,
//                   G4double alpha,
// 		     G4int nLaguerre  ) 
//
//
// -------------------------------------------------------------------------
//
// Gauss-Laguerre method for integration of pow(x,alpha)*exp(-x)*pFunction(x)
// from zero up to infinity. pFunction is evaluated in fNumber points for which
// fAbscissa[i] and fWeight[i] arrays were created in constructor
//
// G4double  Integral() const 

// ------------------------------- HISTORY --------------------------------
//
// 13.05.97 V.Grichine (Vladimir.Grichine@cern.chz0

#ifndef G4GAUSSLAGUERREQ_HH
#define G4GAUSSLAGUERREQ_HH

#include "G4VGaussianQuadrature.hh"

class G4GaussLaguerreQ : public G4VGaussianQuadrature
{
public:
        G4GaussLaguerreQ( function pFunction,
	                  G4double alpha,
			  G4int nLaguerre           ) ;
			       
        // Methods
			     
        G4double Integral() const ;


protected:

private:

} ;


#endif
