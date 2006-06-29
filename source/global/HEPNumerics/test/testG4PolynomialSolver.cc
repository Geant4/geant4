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
// $Id: testG4PolynomialSolver.cc,v 1.6 2006-06-29 19:00:35 gunter Exp $
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
