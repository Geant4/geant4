// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// 
// class G4PolynomialSolver
//
// Implementation
//
// 19.12.00 E.Medernach, First implementation
//

// Class description:
//
//   This class allows the user to solve a polynomial equation
//   with a great precision. This is used by Implicit Equation solver.
//
//   For the moment, we use a Bezier clipping method to solve the polynomial.
//

// How to use it:
//
//   You have to create a class that is the function to be solved.
//   This class could have internal parameters to allow you to change
//   the equation to be solved without recreating a new one.
//
//   You define a Polynomial solver like this one:
//   G4PolynomialSolver<MyFunctionClass,G4double(MyFunctionClass::*)(G4double)>
//     PolySolver (&MyFunction,
//                 &MyFunctionClass::Function,
//                 &MyFunctionClass::Derivative,
//                 precision);
//
//   The precision is relative to the function you solve.
//
//   In MyFunctionClass, you have to provide the function to solve and its derivative:
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
//   Then to have a root inside an interval [IntervalMin,IntervalMax] you have to do the following:
//   MyRoot = PolySolver.solve(IntervalMin,IntervalMax);
//

#include  "globals.hh"
     
#ifndef G4POL_SOLVER_HH
#define G4POL_SOLVER_HH



template <class T, class F>
class G4PolynomialSolver 
{
public:
  
  G4PolynomialSolver(T* typeF,F func,F deriv,G4double precision);  
  ~G4PolynomialSolver();
  

  G4double solve (G4double IntervalMin,G4double IntervalMax);
  
private:
  //General Newton method with Bezier Clipping
  G4double Newton (G4double IntervalMin,G4double IntervalMax);
  
  // Works for polynomial of order less or equal than 4.
  // But could be changed to work for polynomial of any order providing
  // that we find the bezier control points.

  //   This is just one iteration of Bezier Clipping
  int BezierClipping(G4double *IntervalMin,G4double *IntervalMax);


  T* FunctionClass ;
  F Function ;
  F Derivative ;
  
  G4double Precision;
};

#include "G4PolynomialSolver.icc"

#endif 
