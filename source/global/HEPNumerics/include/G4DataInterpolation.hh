// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DataInterpolation.hh,v 1.2 1999-11-16 17:30:56 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class description:
//
// The class consists of some methods for data interpolations and extrapolations.
// The methods based mainly on recommendations given in the book : An introduction to
// NUMERICAL METHODS IN C++, B.H. Flowers, Claredon Press, Oxford, 1995
// 
// ------------------------------ Data members: ---------------------------------
//
//              fArgument and fFunction - pointers to data table to be interpolated
//                                        for y[i] and x[i] respectively
//              fNumber - the corresponding table size
// ......
// G4DataInterpolation( G4double pX[], G4double pY[], G4int number    )
//
// Constructor for initializing of fArgument, fFunction and fNumber data members:
// ......
// G4DataInterpolation( G4double pX[], G4double pY[], G4int number,
//			G4double pFirstDerStart, G4double pFirstDerFinish  ) 
//
// Constructor for cubic spline interpolation. It creates the array 
// fSecondDerivative[0,...fNumber-1] which is used in this interpolation by
// the function: 
// ....
// ~G4DataInterpolation() 
//
// Destructor deletes dynamically created arrays for data members: fArgument,
// fFunction and fSecondDerivative, all have dimension of fNumber
//
// ------------------------------ Methods: ----------------------------------------
//
// G4double PolynomInterpolation(G4double pX, G4double& deltaY ) const
//
// This function returns the value P(pX), where P(x) is polynom of fNumber-1 degree
// such that P(fArgument[i]) = fFunction[i], for i = 0, ..., fNumber-1  .
// ........
// void PolIntCoefficient( G4double cof[]) const 
//
// Given arrays fArgument[0,..,fNumber-1] and fFunction[0,..,fNumber-1] , this
// function calculates an array of coefficients. The coefficients don't provide
// usually (fNumber>10) better accuracy for polynom interpolation, as compared with
// PolynomInterpolation function. They could be used instead for derivate 
// calculations and some other applications.
// .........
// G4double RationalPolInterpolation(G4double pX, G4double& deltaY ) const 
//
// The function returns diagonal rational function (Bulirsch and Stoer algorithm
// of Neville type) Pn(x)/Qm(x) where P and Q are polynoms.
// Tests showed the method is not stable and hasn't advantage if compared with
// polynomial interpolation
// ................
// G4double CubicSplineInterpolation(G4double pX) const 
//
// Cubic spline interpolation in point pX for function given by the table:
// fArgument, fFunction. The constructor, which creates fSecondDerivative, must be
// called before. The function works optimal, if sequential calls are in random
// values of pX.
// ..................
// G4double FastCubicSpline(G4double pX, G4int index) const 
//
// Return cubic spline interpolation in the point pX which is located between
// fArgument[index] and fArgument[index+1]. It is usually called in sequence of
// known from external analysis values of index.
// .........
// G4int LocateArgument(G4double pX) const 
//
// Given argument pX, returns index k, so that pX bracketed by fArgument[k] and
// fArgument[k+1]
// ......................
// void CorrelatedSearch( G4double pX, G4int& index ) const 
//
// Given a value pX, returns a value 'index' such that pX is between fArgument[index]
// and fArgument[index+1]. fArgument MUST BE MONOTONIC, either increasing or
// decreasing. If index = -1 or fNumber, this indicates that pX is out of range.
// The value index on input is taken as the initial approximation for index on
// output.

// --------------------------------- History: --------------------------------------
//
//  3.4.97 V.Grichine (Vladimir.Grichine@cern.ch)
//


#ifndef G4DATAINTERPOLATION_HH
#define G4DATAINTERPOLATION_HH

#include "globals.hh"

class G4DataInterpolation
{
public:
            G4DataInterpolation( G4double pX[], 
	                         G4double pY[], 
			         G4int number    );

// Constructor for cubic spline interpolation. It creates fSecond Deivative array 
// as well as fArgument and fFunction
				 
	    G4DataInterpolation( G4double pX[], 
	                         G4double pY[], 
			         G4int number,
				 G4double pFirstDerStart,
				 G4double pFirstDerFinish  ) ;

           ~G4DataInterpolation() ;
	   
	    G4double PolynomInterpolation( G4double pX,
	                                   G4double& deltaY ) const ;
	    
	    void PolIntCoefficient( G4double cof[]) const ;

            G4double RationalPolInterpolation( G4double pX,
	                                       G4double& deltaY ) const ;

	    G4double CubicSplineInterpolation( G4double pX ) const ;				      

            G4double FastCubicSpline( G4double pX, 
				      G4int index ) const ;

	    G4int LocateArgument( G4double pX ) const ;
	    
	    void CorrelatedSearch( G4double pX,
	                           G4int& index ) const ;
	                   
protected:

private:
	   G4double* fArgument ;
           G4double* fFunction ;
	   G4double* fSecondDerivative ;
           G4int     fNumber ;
} ;

#endif
