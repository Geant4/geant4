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
// $Id: G4PolynomialSolver.hh 67970 2013-03-13 10:10:06Z gcosmo $
// 
// class G4PolynomialSolver
//
// Class description:
//
//   G4PolynomialSolver allows the user to solve a polynomial equation
//   with a great precision. This is used by Implicit Equation solver.
//
//   The Bezier clipping method is used to solve the polynomial.
//
// How to use it:
//   Create a class that is the function to be solved.
//   This class could have internal parameters to allow to change
//   the equation to be solved without recreating a new one.
//
//   Define a Polynomial solver, example:
//   G4PolynomialSolver<MyFunctionClass,G4double(MyFunctionClass::*)(G4double)>
//     PolySolver (&MyFunction,
//                 &MyFunctionClass::Function,
//                 &MyFunctionClass::Derivative,
//                 precision);
//
//   The precision is relative to the function to solve.
//
//   In MyFunctionClass, provide the function to solve and its derivative:
//   Example of function to provide :
//
//   x,y,z,dx,dy,dz,Rmin,Rmax are internal variables of MyFunctionClass
//
//   G4double MyFunctionClass::Function(G4double value)
//   {
//     G4double Lx,Ly,Lz;
//     G4double result;  
//   
//     Lx = x + value*dx;
//     Ly = y + value*dy;
//     Lz = z + value*dz;
//   
//     result = TorusEquation(Lx,Ly,Lz,Rmax,Rmin);
//     
//     return result ;  
//   }    
// 
//   G4double MyFunctionClass::Derivative(G4double value)
//   {
//     G4double Lx,Ly,Lz;
//     G4double result;  
//     
//     Lx = x + value*dx;
//     Ly = y + value*dy;
//     Lz = z + value*dz;
//      
//     result = dx*TorusDerivativeX(Lx,Ly,Lz,Rmax,Rmin);
//     result += dy*TorusDerivativeY(Lx,Ly,Lz,Rmax,Rmin);
//     result += dz*TorusDerivativeZ(Lx,Ly,Lz,Rmax,Rmin);
//   
//     return result;
//   }
//   
//   Then to have a root inside an interval [IntervalMin,IntervalMax] do the
//   following:
//
//   MyRoot = PolySolver.solve(IntervalMin,IntervalMax);
//

// History:
//
// - 19.12.00 E.Medernach, First implementation
//

#ifndef G4POL_SOLVER_HH
#define G4POL_SOLVER_HH

#include  "globals.hh"

template <class T, class F>
class G4PolynomialSolver 
{
public:  // with description
  
  G4PolynomialSolver(T* typeF, F func, F deriv, G4double precision);  
  ~G4PolynomialSolver();
  

  G4double solve (G4double IntervalMin, G4double IntervalMax);
  
private:

  G4double Newton (G4double IntervalMin, G4double IntervalMax);
    //General Newton method with Bezier Clipping

  // Works for polynomial of order less or equal than 4.
  // But could be changed to work for polynomial of any order providing
  // that we find the bezier control points.

  G4int BezierClipping(G4double *IntervalMin, G4double *IntervalMax);
    //   This is just one iteration of Bezier Clipping


  T* FunctionClass ;
  F Function ;
  F Derivative ;
  
  G4double Precision;
};

#include "G4PolynomialSolver.icc"

#endif 
