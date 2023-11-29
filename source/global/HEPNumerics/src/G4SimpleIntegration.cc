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
// G4SimpleIntegration class implementation
//
// Author: V.Grichine, 26.03.1997
// --------------------------------------------------------------------

#include "G4SimpleIntegration.hh"
#include "globals.hh"

G4SimpleIntegration::G4SimpleIntegration(function pFunction)
  : fFunction(pFunction)
{}

G4SimpleIntegration::G4SimpleIntegration(function pFunction,
                                         G4double pTolerance)
  : fFunction(pFunction)
  , fTolerance(pTolerance)
{}

// Simple integration methods

G4double G4SimpleIntegration::Trapezoidal(G4double xInitial, G4double xFinal,
                                          G4int iterationNumber)
{
  G4double Step = (xFinal - xInitial) / iterationNumber;
  G4double mean = (fFunction(xInitial) + fFunction(xFinal)) * 0.5;
  G4double x    = xInitial;
  for(G4int i = 1; i < iterationNumber; ++i)
  {
    x += Step;
    mean += fFunction(x);
  }
  return mean * Step;
}

G4double G4SimpleIntegration::MidPoint(G4double xInitial, G4double xFinal,
                                       G4int iterationNumber)
{
  G4double Step = (xFinal - xInitial) / iterationNumber;
  G4double x    = xInitial + 0.5 * Step;
  G4double mean = fFunction(x);
  for(G4int i = 1; i < iterationNumber; ++i)
  {
    x += Step;
    mean += fFunction(x);
  }
  return mean * Step;
}

G4double G4SimpleIntegration::Gauss(G4double xInitial, G4double xFinal,
                                    G4int iterationNumber)
{
  G4double x                 = 0.;
  static const G4double root = 1.0 / std::sqrt(3.0);
  G4double Step              = (xFinal - xInitial) / (2.0 * iterationNumber);
  G4double delta             = Step * root;
  G4double mean              = 0.0;
  for(G4int i = 0; i < iterationNumber; ++i)
  {
    x = (2 * i + 1) * Step;
    mean += (fFunction(x + delta) + fFunction(x - delta));
  }
  return mean * Step;
}

G4double G4SimpleIntegration::Simpson(G4double xInitial, G4double xFinal,
                                      G4int iterationNumber)
{
  G4double Step  = (xFinal - xInitial) / iterationNumber;
  G4double x     = xInitial;
  G4double xPlus = xInitial + 0.5 * Step;
  G4double mean  = (fFunction(xInitial) + fFunction(xFinal)) * 0.5;
  G4double sum   = fFunction(xPlus);
  for(G4int i = 1; i < iterationNumber; ++i)
  {
    x += Step;
    xPlus += Step;
    mean += fFunction(x);
    sum += fFunction(xPlus);
  }
  mean += 2.0 * sum;
  return mean * Step / 3.0;
}

// Adaptive Gauss integration

G4double G4SimpleIntegration::AdaptGaussIntegration(G4double xInitial,
                                                    G4double xFinal)
{
  G4int depth  = 0;
  G4double sum = 0.0;
  AdaptGauss(xInitial, xFinal, sum, depth);
  return sum;
}

G4double G4SimpleIntegration::Gauss(G4double xInitial, G4double xFinal)
{
  static const G4double root = 1.0 / std::sqrt(3.0);

  G4double xMean = (xInitial + xFinal) / 2.0;
  G4double Step  = (xFinal - xInitial) / 2.0;
  G4double delta = Step * root;
  G4double sum   = (fFunction(xMean + delta) + fFunction(xMean - delta));

  return sum * Step;
}

void G4SimpleIntegration::AdaptGauss(G4double xInitial, G4double xFinal,
                                     G4double& sum, G4int& depth)
{
  if(depth > fMaxDepth)
  {
    G4Exception("G4SimpleIntegration::AdaptGauss()", "Error", FatalException,
                "Function varies too rapidly !");
  }
  G4double xMean     = (xInitial + xFinal) / 2.0;
  G4double leftHalf  = Gauss(xInitial, xMean);
  G4double rightHalf = Gauss(xMean, xFinal);
  G4double full      = Gauss(xInitial, xFinal);
  if(std::fabs(leftHalf + rightHalf - full) < fTolerance)
  {
    sum += full;
  }
  else
  {
    ++depth;
    AdaptGauss(xInitial, xMean, sum, depth);
    AdaptGauss(xMean, xFinal, sum, depth);
  }
}
