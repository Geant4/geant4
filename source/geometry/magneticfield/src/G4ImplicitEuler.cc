// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ImplicitEuler.cc,v 1.3 2000-11-01 15:15:53 gcosmo Exp $
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
G4MagErrorStepper(EqRhs, numberOfVariables),
  fNumberOfVariables(numberOfVariables)
{
}


////////////////////////////////////////////////////////////////////////
//
// Destructor

G4ImplicitEuler::~G4ImplicitEuler()
{
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
  //  const G4int nvar = 6 ;
  G4double* dydxTemp = new G4double[fNumberOfVariables] ;
  G4double* yTemp    = new G4double[fNumberOfVariables] ;

  G4int i;

  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yTemp[i] = yIn[i] + h*dydx[i] ;          
  }
  
  RightHandSide(yTemp,dydxTemp);
  
  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yOut[i] = yIn[i] + 0.5 * h * ( dydx[i] + dydxTemp[i] );
  }

  // NormaliseTangentVector( yOut );           
  
  return ;
}  
