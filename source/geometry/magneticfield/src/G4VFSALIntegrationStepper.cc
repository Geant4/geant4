

#include "G4VFSALIntegrationStepper.hh"

// Constructor for stepper abstract base class. 
// 

G4VFSALIntegrationStepper::G4VFSALIntegrationStepper(G4EquationOfMotion* Equation,
					       G4int       num_integration_vars,
					       G4int       num_state_vars)
  : fEquation_Rhs(Equation),
    fNoIntegrationVariables(num_integration_vars),
    fNoStateVariables(num_state_vars)
    // fNumberOfVariables( std::max(num_var,fNoStateVariables) )
{
}

G4VFSALIntegrationStepper::~G4VFSALIntegrationStepper()
{
}

void G4VFSALIntegrationStepper::ComputeRightHandSide( const G4double y[], G4double dydx[] ) 
{
  this->RightHandSide( y, dydx );
//	fEquation_Rhs->RightHandSide(y, dydx);
//        increasefNORHSCalls();
}

void G4VFSALIntegrationStepper::increasefNORHSCalls(){
    //    std::cout<<"Yeah, I was called!";
    fNoRHSCalls++;
}


void G4VFSALIntegrationStepper::RightHandSide( const  double y[], double dydx[] )
{
    fEquation_Rhs-> RightHandSide(y, dydx);
    increasefNORHSCalls();
}


