// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HelixExplicitEuler.cc,v 1.3 2000-04-12 18:29:26 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4HelixExplicitEuler.hh"
#include "G4ThreeVector.hh"

//
//  Helix Explicit Euler: x_1 = x_0 + helix(h)
//  with helix(h) being a helix piece of length h
//  W.Wander <wwc@mit.edu> 12/09/97 
//

// -------------------------------------------------------------------------

// most simple approach for solving linear differential equations.
// Take the current derivative and add it to the current position.
//

void
G4HelixExplicitEuler::DumbStepper( const G4double  yIn[],
				   G4ThreeVector   Bfld,
				   G4double  h,
				   G4double  yOut[])
{
  AdvanceHelix(yIn, Bfld, h, yOut);

  // NormaliseTangentVector( yOut );  // this could harm more than it helps 
  return ;
}  
