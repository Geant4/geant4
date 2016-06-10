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
// $Id: G4VGaussianQuadrature.cc 69546 2013-05-08 09:50:34Z gcosmo $
//
// Implementation file for G4VGaussianQuadrature virtual base class
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4VGaussianQuadrature.hh"

G4VGaussianQuadrature::G4VGaussianQuadrature( function pFunction )
  : fFunction(pFunction), fAbscissa(0), fWeight(0), fNumber(0)
{
}

// -------------------------------------------------------------------
//
// Virtual destructor which deletes dynamically allocated memory
//

G4VGaussianQuadrature::~G4VGaussianQuadrature() 
{
   delete[] fAbscissa ;
   delete[] fWeight   ;
}

// -------------------------- Access functions ----------------------------------

G4double
G4VGaussianQuadrature::GetAbscissa(G4int index) const
{
   return fAbscissa[index] ;
}

G4double
G4VGaussianQuadrature::GetWeight(G4int index) const
{
   return fWeight[index] ;
}

G4int G4VGaussianQuadrature::GetNumber() const
{
   return fNumber ;
}

// ----------------------------------------------------------------------------
//
// Auxiliary function which returns the value of std::log(gamma-function(x))
//

G4double 
G4VGaussianQuadrature::GammaLogarithm(G4double xx)
{

// Returns the value ln(Gamma(xx) for xx > 0.  Full accuracy is obtained for 
// xx > 1. For 0 < xx < 1. the reflection formula (6.1.4) can be used first.
// (Adapted from Numerical Recipes in C)

  static const G4double cof[6] = { 76.18009172947146,     -86.50532032941677,
                                  24.01409824083091,      -1.231739572450155,
                                   0.1208650973866179e-2, -0.5395239384953e-5 };
  G4double x = xx - 1.0;
  G4double tmp = x + 5.5;
  tmp -= (x + 0.5) * std::log(tmp);
  G4double ser = 1.000000000190015;

  for ( size_t j = 0; j <= 5; j++ )
  {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp + std::log(2.5066282746310005*ser);
}
