// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGaussianQuadrature.cc,v 1.2 1999-11-16 17:31:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation file for G4VGaussianQuadrature virtual base class
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4VGaussianQuadrature.hh"




G4VGaussianQuadrature::G4VGaussianQuadrature( function pFunction )
{
   fFunction = pFunction ;
   fAbscissa = 0;
   fWeight = 0;
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


// ----------------------------------------------------------------------------
//
// Auxiliary function which returns the value of log(gamma-function(x))
//

G4double 
G4VGaussianQuadrature::GammaLogarithm(G4double xx)
{

// Returns the value ln(Gamma(xx) for xx > 0.  Full accuracy is obtained for 
// xx > 1. For 0 < xx < 1. the reflection formula (6.1.4) can be used first.
// (Adapted from Numerical Recipes in C)

  static G4double cof[6] = { 76.18009172947146,     -86.50532032941677,
                             24.01409824083091,      -1.231739572450155,
                              0.1208650973866179e-2, -0.5395239384953e-5  } ;
  register HepInt j;
  G4double x = xx - 1.0;
  G4double tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  G4double ser = 1.000000000190015;

  for ( j = 0; j <= 5; j++ )
  {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp + log(2.5066282746310005*ser);
}


   
