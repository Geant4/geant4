// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MagIntegratorStepper.cc,v 2.0 1998/07/02 16:56:24 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4MagIntegratorStepper.hh"

// Constructor for creation of coefficient in Lorentz motion equation 
// as well as initialisation of geometry constants.
// particleCharge is the particle charge in the elementary charge
// momentumXc is the particle momentum multiplied by the speed of light in MeV,
// 

G4MagIntegratorStepper::G4MagIntegratorStepper(G4Mag_EqRhs *EqRhs)
{
   theEquation_Rhs= EqRhs;
}

