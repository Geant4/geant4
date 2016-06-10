//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4SimpleRunge.cc 66356 2012-12-18 09:02:32Z gcosmo $
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
