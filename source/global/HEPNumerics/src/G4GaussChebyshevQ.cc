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
// G4GaussChebyshevQ class implementation
//
// Author: V.Grichine, 13.05.1997
// --------------------------------------------------------------------

#include "G4GaussChebyshevQ.hh"
#include "G4PhysicalConstants.hh"

// -----------------------------------------------------
//
// Constructor for Gauss-Chebyshev quadrature method
//
G4GaussChebyshevQ::G4GaussChebyshevQ(function pFunction, G4int nChebyshev)
  : G4VGaussianQuadrature(pFunction)
{
  fNumber      = nChebyshev;  // Try to reduce fNumber twice ??
  G4double cof = pi / fNumber;
  fAbscissa    = new G4double[fNumber];
  fWeight      = new G4double[fNumber];
  for(G4int i = 0; i < fNumber; ++i)
  {
    fAbscissa[i] = std::cos(cof * (i + 0.5));
    fWeight[i]   = cof * std::sqrt(1 - fAbscissa[i] * fAbscissa[i]);
  }
}

// -------------------------------------------------------------------------------
//
// Integrates function pointed by fFunction from a to b by Gauss-Chebyshev
// quadrature method
//
G4double G4GaussChebyshevQ::Integral(G4double a, G4double b) const
{
  G4double xDiff = 0.5 * (b - a), xMean = 0.5 * (a + b), dx = 0.0,
           integral = 0.0;

  for(G4int i = 0; i < fNumber; ++i)
  {
    dx = xDiff * fAbscissa[i];
    integral += fWeight[i] * fFunction(xMean + dx);
  }
  return integral *= xDiff;
}
