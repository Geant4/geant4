// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixSimpleRunge.cc,v 2.3 1998/11/13 14:30:22 japost Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4HelixSimpleRunge.hh"
#include "G4ThreeVector.hh"

//
//  Simple Runge:
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//  W.Wander <wwc@mit.edu> 12/09/97 
//

// -------------------------------------------------------------------------

// second order solver
// take the derivative at a position to be assumed at the middle of the
// Step and add it to the current position.
//

void
G4HelixSimpleRunge::DumbStepper( const G4double  yIn[],
				 const G4double  dydx[],
				 const G4double  h,
				 G4double  yOut[])
{
  const G4int nvar = 6 ;
  G4double dydxTemp[nvar];
  G4double yTemp[nvar];   // , yAdd[nvar];

  AdvanceHelix( yIn, dydx, 0.5 * h, yTemp);
  
  RightHandSide(yTemp,dydxTemp);

  AdvanceHelix( yIn, dydxTemp, h, yOut);
  
  // NormaliseTangentVector( yOut );           
  
  return ;
}  
