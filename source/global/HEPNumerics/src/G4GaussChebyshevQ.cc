//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GaussChebyshevQ.cc,v 1.3 2001-07-11 10:00:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4GaussChebyshevQ.hh"

// -----------------------------------------------------
//
// Constructor for Gauss-Chebyshev quadrature method

G4GaussChebyshevQ::G4GaussChebyshevQ( function pFunction ,
			              G4int nChebyshev       )
   : G4VGaussianQuadrature(pFunction)
{
   G4int i ;
   fNumber = nChebyshev  ;   // Try to reduce fNumber twice ??
   G4double cof = pi/fNumber ;
   fAbscissa = new G4double[fNumber] ;
   fWeight = new G4double[fNumber] ;
   for(i=0;i<fNumber;i++)
   {
      fAbscissa[i] = cos(cof*(i + 0.5)) ;
      fWeight[i] = cof*sqrt(1 - fAbscissa[i]*fAbscissa[i]) ;
   }
}

// ----------------------------------------------------------------------
//

G4GaussChebyshevQ::~G4GaussChebyshevQ()
{
   ;
}

// -------------------------------------------------------------------------------
//
// Integrates function pointed by fFunction from a to b by Gauss-Chebyshev 
// quadrature method

G4double 
G4GaussChebyshevQ::Integral(G4double a, G4double b) const 
{
   G4int i ;
   G4double xDiff, xMean, dx, integral = 0.0 ;
   
   xMean = 0.5*(a + b) ;
   xDiff = 0.5*(b - a) ;
   for(i=0;i<fNumber;i++)
   {
      dx = xDiff*fAbscissa[i] ;
      integral += fWeight[i]*fFunction(xMean + dx)  ;
   }
   return integral *= xDiff ;
}
