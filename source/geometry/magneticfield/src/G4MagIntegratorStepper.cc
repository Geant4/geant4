// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagIntegratorStepper.cc,v 1.4 2000-11-01 15:15:53 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4MagIntegratorStepper.hh"

// Constructor for stepper abstract base class. 
// 

G4MagIntegratorStepper::G4MagIntegratorStepper(G4Mag_EqRhs *EqRhs,
					       G4int       num_var)
  : fEquation_Rhs(EqRhs),
    fNumberOfVariables(num_var)
{
}

G4MagIntegratorStepper::~G4MagIntegratorStepper()
{
}

void G4MagIntegratorStepper::RightHandSide( const  double y[], double dydx[] )   
{
  fEquation_Rhs-> RightHandSide(y, dydx);
}
