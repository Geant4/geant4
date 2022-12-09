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
// G4SimpleIntegration
//
// Class description:
//
// Class for realisation of simple numerical methodes for integration of
// functions with signature: double f(double). The methods based mainly on
// algorithms given in the book:
//   An introduction to NUMERICAL METHODS IN C++,
//   B.H. Flowers, Claredon Press, Oxford, 1995.

// Author: V.Grichine, 26.03.1997
// --------------------------------------------------------------------
#ifndef G4SIMPLEINTEGRATION_HH
#define G4SIMPLEINTEGRATION_HH 1

#include "G4Types.hh"

using function = G4double (*)(G4double);

class G4SimpleIntegration
{
 public:
  explicit G4SimpleIntegration(function pFunction);

  G4SimpleIntegration(function pFunction, G4double pTolerance);

  ~G4SimpleIntegration() = default;

  G4SimpleIntegration(const G4SimpleIntegration&) = delete;
  G4SimpleIntegration& operator=(const G4SimpleIntegration&) = delete;
  // Private copy constructor and assignment operator.

  // Simple integration methods:
  // Trapezoidal, MidPoint, Gauss and Simpson(double a,double b,int n)
  // - integrate function pointed by fFunction from a to b by n iterations,
  //   i.e. with Step (b-a)/n according to the correspondent method.

  G4double Trapezoidal(G4double xInitial, G4double xFinal,
                       G4int iterationNumber);

  G4double MidPoint(G4double xInitial, G4double xFinal, G4int iterationNumber);

  G4double Gauss(G4double xInitial, G4double xFinal, G4int iterationNumber);

  G4double Simpson(G4double xInitial, G4double xFinal, G4int iterationNumber);

  // Adaptive Gauss integration with accuracy ~ fTolerance

  G4double AdaptGaussIntegration(G4double xInitial, G4double xFinal);
  // Integrate function from a to be with accuracy <= fTolerance

 protected:
  G4double Gauss(G4double xInitial, G4double xFinal);

  void AdaptGauss(G4double xInitial, G4double xFinal, G4double& sum,
                  G4int& depth);

 private:
  function fFunction;              // pointer to the function to be integrated
  G4double fTolerance   = 0.0001;  // accuracy of integration
  const G4int fMaxDepth = 100;     // constant maximum iteration depth
};

#endif
