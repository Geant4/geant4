//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ImplicitEuler.cc,v 1.5 2002-11-29 13:45:21 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Implicit Euler:
//
//        x_1 = x_0 + h/2 * ( dx(t_0,x_0) + dx(t_0+h,x_0+h*dx(t_0,x_0) ) )
//
// second order solver
// Take the current derivative and add it to the current position.
// Take the output and its derivative. Add the mean of both derivatives
// to form the final output
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// 6.11.98 V.Grichine, new constructor, fNumberOfVariables
//

#include "G4ImplicitEuler.hh"
#include "G4ThreeVector.hh"

/////////////////////////////////////////////////////////////////////////
//
// Constructor

G4ImplicitEuler::G4ImplicitEuler(G4Mag_EqRhs *EqRhs, 
                                 G4int numberOfVariables): 
G4MagErrorStepper(EqRhs, numberOfVariables)
{
  unsigned int noVariables= G4std::max(numberOfVariables,8); // For Time .. 7+1
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
