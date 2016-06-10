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
// $Id: G4ChebyshevApproximation.hh 67970 2013-03-13 10:10:06Z gcosmo $
//
// Class description:
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
//                           G4double b       ) 
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
//                           G4int m,
//                           G4double a,
//                           G4double b       ) 
//
// ------------------------------------------------------
//
// Constructor for creation of Chebyshev coefficients for integral
// from pFunction.
//
// G4ChebyshevApproximation( function pFunction,
//                           G4double a,
//                           G4double b, 
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

// --------------------------- HISTORY --------------------------------------
//
//  24.04.97   V.Grichine ( Vladimir.Grichine@cern.ch )

#ifndef G4CHEBYSHEVAPPROXIMATION_HH
#define G4CHEBYSHEVAPPROXIMATION_HH

#include "globals.hh"

typedef G4double (*function)(G4double) ;

class G4ChebyshevApproximation
{
  public:  // with description

    G4ChebyshevApproximation( function pFunction,
                              G4int n, 
                              G4double a,
                              G4double b       ) ;
      //
      // Constructor for creation of Chebyshev coefficients for m-derivative
      // from pFunction. The value of m ! MUST BE ! < n , because the result
      // array of fChebyshevCof will be of (n-m) size.

    G4ChebyshevApproximation( function pFunction,
                              G4int n,
                              G4int m,
                              G4double a,
                              G4double b       ) ;
      //
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
       
  private:

    G4ChebyshevApproximation(const G4ChebyshevApproximation&);
    G4ChebyshevApproximation& operator=(const G4ChebyshevApproximation&);

  private:

    function   fFunction ;
    G4int      fNumber ;
    G4double*  fChebyshevCof ;
    G4double   fMean ;
    G4double   fDiff ;
};

#endif
