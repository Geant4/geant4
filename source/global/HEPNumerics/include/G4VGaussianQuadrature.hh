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
// G4VGaussianQuadrature
//
// Class description:
//
// Base Class for realisation of numerical methodes for integration of functions
// with signature double f(double) by Gaussian quadrature methods
// Roots of ortogonal polynoms and corresponding weights are calculated based on
// iteration method (by bisection Newton algorithm). Constant values for initial
// approximations were derived from the book:
//   M. Abramowitz, I. Stegun, Handbook of mathematical functions,
//   DOVER Publications INC, New York 1965 ; chapters 9, 10, and 22.

// Author: V.Grichine, 18.04.1997
// --------------------------------------------------------------------
#ifndef G4VGAUSSIANQUADRATURE_HH
#define G4VGAUSSIANQUADRATURE_HH 1

#include "globals.hh"

using function = G4double (*)(G4double);

class G4VGaussianQuadrature
{
 public:
  explicit G4VGaussianQuadrature(function pFunction);
  // Base constructor

  virtual ~G4VGaussianQuadrature();
  // Virtual destructor

  G4VGaussianQuadrature(const G4VGaussianQuadrature&) = delete;
  G4VGaussianQuadrature& operator=(const G4VGaussianQuadrature&) = delete;

  G4double GetAbscissa(G4int index) const;
  G4double GetWeight(G4int index) const;
  G4int GetNumber() const;
  // Access functions

 protected:
  G4double GammaLogarithm(G4double xx);
  // Auxiliary function which returns the value of std::log(gamma-function(x))

  //  Data members common for GaussianQuadrature family
  //
  function fFunction;             // pointer to the function to be integrated
  G4double* fAbscissa = nullptr;  // array of abscissas
  G4double* fWeight   = nullptr;  // array of corresponding weights
  G4int fNumber = 0;  // the number of points in fAbscissa and fWeight arrays
};

#endif
