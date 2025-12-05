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
// G4FunctionSolver
//
// Class description:
//
// Templated utility class to solve equation F(x) = 0 on interval [a, b]
// Parameters of the class: tolerance - relative accuracy of the method,
// maxIter - max number of iterations. User should provide a class with
// only one method: T_Function::Function(G4double x). The main method
// of the class is 'G4bool G4FunctionSolver::FindRoot(G4double& x)', where
// user defines initial value x, which is modifided to be a final solution
// of the equation.
//
// If the equation cannot be resolved with required accuracy FindRoot()
// method returns "false".
//
// Parameters of the solver may be changed before the new call via Set
// methods.
//
// Created 04.10.2025 V.Ivanchenko
//  

#ifndef G4FunctionSolver_h
#define G4FunctionSolver_h 1

#include "globals.hh"

template <class T_Function> class G4FunctionSolver 
{
public:
	
  G4FunctionSolver(T_Function* ff, const G4int iterations, const G4double tol)
    : maxIter(iterations), tolerance(tol), tF(ff) {};

  // copy constructor	
  G4FunctionSolver(const G4FunctionSolver& right) = delete;

  // destructor
  ~G4FunctionSolver() = default;
	
  // operators
  G4FunctionSolver& operator=(const G4FunctionSolver& right) = delete;
  G4bool operator==(const G4FunctionSolver& right) const = delete;
  G4bool operator!=(const G4FunctionSolver& right) const = delete;

  inline void SetMaxIterations(const G4int iterations)
  { maxIter = iterations; }

  inline void SetTolerance(const G4double epsilon)
  { tolerance = epsilon; }

  inline void SetIntervalLimits(const G4double Limit1, const G4double Limit2)
  {
    aa = std::min(Limit1, Limit2);
    bb = std::max(Limit1, Limit2);
  }

  // Calculates the root of the equation Function(x)=0
  inline G4bool FindRoot(G4double& x)
  {
    G4double a = aa;
    G4double b = bb;

    // check the interval before the start
    x = std::min(std::max(x, a), b);

    // check initial function
    G4double fc = tF->Function(x);
    if (0.0 == fc) { return true; }

    // define accuracy in X
    G4double epsX = tolerance*x;
    // define accuracy in Y
    G4double epsY = tolerance*fc;

    // the interval is too small
    if (std::abs(a - b) <= epsX)
    {
      x = 0.5*(a + b);
      return true;
    }

    // check edges
    G4double fa = tF->Function(a);
    G4double fb = tF->Function(b);

    // root should be inside interval
    if (fa*fb >= 0.0)
    {
      x = (std::abs(fa) <= std::abs(fb)) ? a : b;
      return (std::min(std::abs(fa), std::abs(fb)) < epsY);
    }

    // fa*fb < 0.0 - finding the root by iterative procedure, 
    // the loop is completed if function is below epsY
    // or if x become close to edges with accuracy epsX 
    for (G4int i = 0; i < maxIter; ++i)
    {
      x = (a*fb - b*fa)/(fb - fa);
      fc = tF->Function(x);
      if (std::abs(fc) < epsY) { return true; }

      G4double delta = std::min((x - a), (b - x));
      if (delta < epsX) { return true; }
      else if (fa*fc < 0.0) { b = x; fb = fc; }
      else { a = x; fa = fc; }
    }
    // number of iterations exceed the limit
    return false;
  }

private:

  // maximum number of iterations
  G4int maxIter;
  // relative accuracy in X and Y
  G4double tolerance;

  // interval limits [a,b] 
  G4double aa{0.0};
  G4double bb{0.0};

  T_Function* tF;
};

#endif
