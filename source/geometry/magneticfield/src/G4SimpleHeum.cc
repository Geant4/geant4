// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SimpleHeum.cc,v 1.3 1999-12-15 14:49:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  Simple Heum:
//        x_1 = x_0 + h *
//                1/4 * dx(t0,x0)  +
//                3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
//
// third order solver
//
//  W.Wander <wwc@mit.edu> 12/09/97 
//  6.11.98. V.Grichine new data member fNumberOfVariables


#include "G4SimpleHeum.hh"
#include "G4ThreeVector.hh"

///////////////////////////////////////////////////////////////////////////
//
// Constructor

G4SimpleHeum::G4SimpleHeum(G4Mag_EqRhs *EqRhs, G4int num_variables): 
  G4MagErrorStepper(EqRhs, num_variables),
  fNumberOfVariables(num_variables)
{
  dydxTemp  = new G4double[fNumberOfVariables] ; 
  dydxTemp2 = new G4double[fNumberOfVariables] ;
  yTemp     = new G4double[fNumberOfVariables] ; 
  yTemp2    = new G4double[fNumberOfVariables] ;
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4SimpleHeum::~G4SimpleHeum()
{
  delete[] dydxTemp;
  delete[] dydxTemp2;
  delete[] yTemp;
  delete[] yTemp2;
}


//////////////////////////////////////////////////////////////////////
//
//

void
G4SimpleHeum::DumbStepper( const G4double  yIn[],
			    const G4double  dydx[],
			    const G4double  h,
			 	  G4double  yOut[])
{
  //  const G4int nvar = 6 ;

  G4int i;

  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yTemp[i] = yIn[i] + (1.0/3.0) * h *  dydx[i] ;
  }
  
  RightHandSide(yTemp,dydxTemp);

  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yTemp2[i] = yIn[i] + (2.0/3.0) * h * dydxTemp[i] ;
  }

  RightHandSide(yTemp2,dydxTemp2);

  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yOut[i] = yIn[i] + h * ( 0.25 * dydx[i] +
			     0.75 * dydxTemp2[i]);
  }
      
  // NormaliseTangentVector( yOut );           
  
  return ;
}  
