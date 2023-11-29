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
// G4ChebyshevApproximation
//
// Class description:
//
// Class creating the Chebyshev approximation for a function pointed by
// fFunction data member. The Chebyshev polinom approximation provides an
// efficient evaluation of minimax polynomial, which (among all polynomials of
// the same degree) has the smallest maximum deviation from the true function.
// The methods based mainly on recommendations given in the book : An
// introduction to NUMERICAL METHODS IN C++, B.H. Flowers, Claredon Press,
// Oxford, 1995

// Author: V.Grichine, 24.04.1997
// --------------------------------------------------------------------
#ifndef G4CHEBYSHEVAPPROXIMATION_HH
#define G4CHEBYSHEVAPPROXIMATION_HH 1

#include "globals.hh"

using function = G4double (*)(G4double);

class G4ChebyshevApproximation
{
 public:
  G4ChebyshevApproximation(function pFunction, G4int n, G4double a, G4double b);
  // Constructor for creation of Chebyshev coefficients for m-derivative
  // from pFunction. The value of m ! MUST BE ! < n , because the result
  // array of fChebyshevCof will be of (n-m) size.
  // It creates the array fChebyshevCof[0,...,fNumber-1], fNumber = n ;
  // which consists of Chebyshev coefficients describing the function pointed
  // by pFunction. The values a and b fixe the interval of validity of
  // Chebyshev approximation.

  G4ChebyshevApproximation(function pFunction, G4int n, G4int m, G4double a,
                           G4double b);
  // Constructor for creation of Chebyshev coefficients for m-derivative
  // from pFunction. The value of m ! MUST BE ! < n , because the result
  // array of fChebyshevCof will be of (n-m) size. There is a definite
  // dependence between the proper selection of n, m, a and b values to get
  // better accuracy of the derivative value.

  G4ChebyshevApproximation(function pFunction, G4double a, G4double b, G4int n);
  // Constructor for creation of Chebyshev coefficients for integral
  // from pFunction.

  ~G4ChebyshevApproximation();
  // Destructor deletes the array of Chebyshev coefficients

  G4ChebyshevApproximation(const G4ChebyshevApproximation&) = delete;
  G4ChebyshevApproximation& operator=(const G4ChebyshevApproximation&) = delete;
  // Copy constructor and assignment operator not allowed.

  G4double GetChebyshevCof(G4int number) const;
  // Access function for Chebyshev coefficients

  G4double ChebyshevEvaluation(G4double x) const;
  // Evaluate the value of fFunction at the point x via the Chebyshev
  // coefficients fChebyshevCof[0,...,fNumber-1]

  void DerivativeChebyshevCof(G4double derCof[]) const;
  // Returns the array derCof[0,...,fNumber-2], the Chebyshev coefficients
  // of the derivative of the function whose coefficients are fChebyshevCof

  void IntegralChebyshevCof(G4double integralCof[]) const;
  // This function produces the array integralCof[0,...,fNumber-1] , the
  // Chebyshev coefficients of the integral of the function whose coefficients
  // are fChebyshevCof. The constant of integration is set so that the integral
  // vanishes at the point (fMean - fDiff)

 private:
  function fFunction;       // pointer to a function considered
  G4int fNumber;            // number of Chebyshev coefficients
  G4double* fChebyshevCof;  // array of Chebyshev coefficients
  G4double fMean;           // (a+b)/2 - mean point of interval
  G4double fDiff;           // (b-a)/2 - half of the interval value
};

#endif
