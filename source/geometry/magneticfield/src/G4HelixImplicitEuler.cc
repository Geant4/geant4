// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixImplicitEuler.cc,v 1.1 1999-01-07 16:07:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4HelixImplicitEuler.hh"
#include "G4ThreeVector.hh"

//
//  Helix Implicit Euler:
//        x_1 = x_0 + 1/2 * ( helix(h,t_0,x_0)
//                          + helix(h,t_0+h,x_0+helix(h,t0,x0) ) )
//  W.Wander <wwc@mit.edu> 12/09/97 
//

// -------------------------------------------------------------------------

// second order solver
// Take the current derivative and add it to the current position.
// Take the output and its derivative. Add the mean of both derivatives
// to form the final output
//

void
G4HelixImplicitEuler::DumbStepper( const G4double  yIn[],
				   const G4double  dydx[],
				   const G4double  h,
				   G4double  yOut[])
{
  const G4int nvar = 6 ;
  G4double dydxTemp[6];
  G4double yTemp[6], yTemp2[6];

  G4int i;

  // Step forward like in the explicit euler case
  AdvanceHelix( yIn, dydx, h, yTemp);

  // now obtain the new field value at the new point
  RightHandSide(yTemp,dydxTemp);

  // and also advance along a helix for this field value
  AdvanceHelix( yIn, dydxTemp, h, yTemp2);

  // we take the average 
  for( i = 0; i < nvar; i++ ) 
    yOut[i] = 0.5 * ( yTemp[i] + yTemp2[i] );

  // NormaliseTangentVector( yOut );           
  
  return ;
}  
