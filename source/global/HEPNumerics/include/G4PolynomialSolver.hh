// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PolynomialSolver.hh,v 1.2 2001-01-29 09:49:54 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
