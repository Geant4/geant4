//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SimpleRunge.cc,v 1.9 2003/11/05 16:31:49 japost Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//  Simple Runge:
//
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//
// Second order solver.
// Takes the derivative at a position to be assumed at the middle of the
// Step and adds it to the current position.
//
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// -------------------------------------------------------------------

#include "G4SimpleRunge.hh"
#include "G4ThreeVector.hh"

////////////////////////////////////////////////////////////////
//
// Constructor

G4SimpleRunge::G4SimpleRunge(G4EquationOfMotion* EqRhs, G4int numberOfVariables)
  : G4MagErrorStepper(EqRhs, numberOfVariables),
    fNumberOfVariables(numberOfVariables)
{
   
   unsigned int noVariables= std::max(numberOfVariables,
					GetNumberOfStateVariables()); 
                                             // To deal with Time >= 7+1 
   dydxTemp = new G4double[noVariables] ;
   yTemp    = new G4double[noVariables] ;
}


/////////////////////////////////////////////////////////////////
//
// Destructor

G4SimpleRunge::~G4SimpleRunge()
{
   delete[] dydxTemp;
   delete[] yTemp;
}

//////////////////////////////////////////////////////////////////
//
//

void
G4SimpleRunge::DumbStepper( const G4double  yIn[],
			    const G4double  dydx[],
			          G4double  h,
			 	  G4double  yOut[])
{
  // Initialise time to t0, needed when it is not updated by the integration.
  yTemp[7] = yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

  G4int i;

  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yTemp[i] = yIn[i] + 0.5 * h*dydx[i] ;
  }
  
  RightHandSide(yTemp,dydxTemp);
  
  for( i = 0; i < fNumberOfVariables; i++ ) 
  {
    yOut[i] = yIn[i] + h * ( dydxTemp[i] );
  }
}  
