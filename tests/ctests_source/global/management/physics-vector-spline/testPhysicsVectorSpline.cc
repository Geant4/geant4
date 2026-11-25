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
// Regression test for G4PhysicsVector spline interpolation.

#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLinearVector.hh"
#include "G4PhysicsLogVector.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
G4double Polynomial(G4double x)
{
  return 1.0 + x * (2.0 + x * (3.0 + x));
}

G4double Derivative(G4double x)
{
  return 2.0 + x * (6.0 + 3.0 * x);
}

G4bool Check(G4PhysicsVector& vec, G4SplineType type)
{
  const auto n = vec.GetVectorLength();
  for (std::size_t i = 0; i < n; ++i)
  {
    vec.PutValue(i, Polynomial(vec.Energy(i)));
  }
  vec.FillSecondDerivatives(type, Derivative(vec.Energy(0)), Derivative(vec.Energy(n - 1)));

  G4bool pass = true;
  for (std::size_t i = 0; i + 1 < n; ++i)
  {
    for (G4double fraction : {0.25, 0.5, 0.75})
    {
      const G4double x = vec.Energy(i) + fraction * (vec.Energy(i + 1) - vec.Energy(i));
      const G4double expected = Polynomial(x);
      const G4double actual = vec.Value(x);
      if (!std::isfinite(actual)
          || std::abs(actual - expected) > 1.e-10 * std::max(1.0, std::abs(expected)))
      {
        std::cerr << "Spline type " << static_cast<G4int>(type) << ", x = " << x << ": expected "
                  << expected << ", got " << actual << '\n';
        pass = false;
      }
    }
  }
  return pass;
}
}  // namespace

int main()
{
  G4bool pass = true;
  for (auto type : {G4SplineType::Base, G4SplineType::FixedEdges})
  {
    G4PhysicsLinearVector linear(1.0, 6.0, 5, true);
    G4PhysicsLogVector log(1.0, 16.0, 4, true);
    const std::vector<G4double> energies{1.0, 1.5, 3.0, 4.0, 7.0, 9.0};
    G4PhysicsFreeVector free(energies, std::vector<G4double>(energies.size()), true);
    pass = Check(linear, type) && pass;
    pass = Check(log, type) && pass;
    pass = Check(free, type) && pass;
  }

  // The clamped spline also supports four nodes.
  G4PhysicsLinearVector four(1.0, 4.0, 3, true);
  pass = Check(four, G4SplineType::FixedEdges) && pass;
  return pass ? 0 : 1;
}
