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
// $Id: G4ImplicitEuler.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
//
//  Implicit Euler:
//
//        x_1 = x_0 + h/2 * ( dx(t_0,x_0) + dx(t_0+h,x_0+h*dx(t_0,x_0) ) )
//
// Second order solver.
// Take the current derivative and add it to the current position.
// Take the output and its derivative. Add the mean of both derivatives
// to form the final output.
//
// W.Wander <wwc@mit.edu> 12/09/97
//
// --------------------------------------------------------------------

#include "G4ImplicitEuler.hh"
#include "G4ThreeVector.hh"

/////////////////////////////////////////////////////////////////////////
//
// Constructor

G4ImplicitEuler::G4ImplicitEuler(G4EquationOfMotion *EqRhs, 
                                 G4int numberOfVariables): 
G4MagErrorStepper(EqRhs, numberOfVariables)
{
  unsigned int noVariables= std::max(numberOfVariables,8); // For Time .. 7+1
  dydxTemp = new G4double[noVariables] ;
  yTemp    = new G4double[noVariables] ;
}


////////////////////////////////////////////////////////////////////////
//
// Destructor

G4ImplicitEuler::~G4ImplicitEuler()
{
  delete[] dydxTemp;
  delete[] yTemp;
}

//////////////////////////////////////////////////////////////////////
//
//

void
G4ImplicitEuler::DumbStepper( const G4double  yIn[],
			      const G4double  dydx[],
			            G4double  h,
			 	    G4double  yOut[])
{
  G4int i;
  const G4int numberOfVariables= GetNumberOfVariables();

  // Initialise time to t0, needed when it is not updated by the integration.
  yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

  for( i = 0; i < numberOfVariables; i++ ) 
  {
    yTemp[i] = yIn[i] + h*dydx[i] ;          
  }
  
  RightHandSide(yTemp,dydxTemp);
  
  for( i = 0; i < numberOfVariables; i++ ) 
  {
    yOut[i] = yIn[i] + 0.5 * h * ( dydx[i] + dydxTemp[i] );
  }

  return ;
}  
