// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExplicitEuler.cc,v 1.2 1999-12-15 14:49:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Explicit Euler: x_1 = x_0 + h * dx_0
//
//  most simple approach for solving linear differential equations.
//  Take the current derivative and add it to the current position.
//
//  W.Wander <wwc@mit.edu> 12/09/97 

#include "G4ExplicitEuler.hh"
#include "G4ThreeVector.hh"

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4ExplicitEuler::G4ExplicitEuler(G4Mag_EqRhs *EqRhs, 
                                 G4int numberOfVariables) : 
                 G4MagErrorStepper(EqRhs, numberOfVariables),
  fNumberOfVariables(numberOfVariables)
{
}


///////////////////////////////////////////////////////////////////////
//
// Destructor

G4ExplicitEuler::~G4ExplicitEuler()
{
   ;
}


///////////////////////////////////////////////////////////////////////
//
//

void
G4ExplicitEuler::DumbStepper( const G4double  yIn[],
			      const G4double  dydx[],
			      const G4double  h,
			 	    G4double  yOut[]        )
{
  //  const G4int nvar = 6 ; 

  G4int i;

  for(i=0;i<fNumberOfVariables;i++)
  {
    yOut[i] = yIn[i] + h*dydx[i] ;             // 1st and only Step 
  }
  // NormaliseTangentVector( yOut );           // this could harm more than
                                               // it helps - FIXME ???
  
  return ;

}  
