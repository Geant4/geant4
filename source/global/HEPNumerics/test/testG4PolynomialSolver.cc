// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4PolynomialSolver.cc,v 1.4 2001-01-29 10:06:12 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4PolynomialSolver.hh"

class MyFunctionClass 
{
public:
  MyFunctionClass (G4double param)
  {
    Parameter = param;
  }
  
  ~MyFunctionClass () {;}

  void setParam(G4double param)
  {
    Parameter = param;    
  }
  

  G4double Function(G4double value)
  {
    G4double result;

    result = value * value - Parameter;
  
    return result ;  
  }    

  G4double Derivative(G4double value)
  {
    G4double result;

    result = 2 * value ;
    
    return result;
  }

private:

  G4double Parameter;

} ;


int main (void)
{
  MyFunctionClass MyFunction(1.0);

  G4PolynomialSolver<MyFunctionClass,G4double(MyFunctionClass::*)(G4double)>
    PolySolver(&MyFunction,
	       &MyFunctionClass::Function,
	       &MyFunctionClass::Derivative,
	       1e-6) ;


  G4double val;

  val = PolySolver.solve(0.0,2.0);

  G4cout << "val = " << val << G4endl;

  MyFunction.setParam(3.0);
  
  val = PolySolver.solve(0.0,2.0);

  G4cout << "val = " << val << G4endl;

  
  return 0;  
}
