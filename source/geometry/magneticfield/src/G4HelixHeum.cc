// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixHeum.cc,v 1.2 1999-12-15 14:49:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4HelixHeum.hh"
#include "G4ThreeVector.hh"

//
//  Simple Heum:
//        x_1 = x_0 + h *
//                1/4 * dx(t0,x0)  +
//                3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
//  W.Wander <wwc@mit.edu> 12/09/97 
//
// -------------------------------------------------------------------------

// third order solver
//

void
G4HelixHeum::DumbStepper( const G4double  yIn[],
			    const G4double  dydx[],
			    const G4double  h,
			 	  G4double  yOut[])
{
  const G4int nvar = 6 ;
  G4double dydxTemp[6], dydxTemp2[6];
  G4double yTemp[6], yAdd1[6], yAdd2[6] , yTemp2[6];

  G4int i;

  AdvanceHelix( yIn, dydx, h, yAdd1 );
  
  AdvanceHelix( yIn, dydx, h/3.0, yTemp );
  RightHandSide(yTemp,dydxTemp);

  AdvanceHelix( yIn, dydxTemp, (2.0 / 3.0) * h, yTemp2 );
  
  RightHandSide(yTemp2,dydxTemp2);

  AdvanceHelix( yIn, dydxTemp2, h, yAdd2 );

  for( i = 0; i < nvar; i++ ) {
    yOut[i] = ( 0.25 * yAdd1[i] + 0.75 * yAdd2[i]);
  }
      
  // NormaliseTangentVector( yOut );           
  
  return ;
}  
