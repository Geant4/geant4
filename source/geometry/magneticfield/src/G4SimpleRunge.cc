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
// $Id: G4SimpleRunge.cc,v 1.5 2002-11-29 13:44:01 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  Simple Runge:
//
//        x_1 = x_0 + h * ( dx( t_0+h/2, x_0 + h/2 * dx( t_0, x_0) ) )
//
// second order solver
// take the derivative at a position to be assumed at the middle of the
// Step and add it to the current position.
//
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// 6.11.98 V.Grichine new data member fNumberOfVariables
//


#include "G4SimpleRunge.hh"
#include "G4ThreeVector.hh"

////////////////////////////////////////////////////////////////
//
// Constructor

G4SimpleRunge::G4SimpleRunge(G4Mag_EqRhs *EqRhs, G4int numberOfVariables)
  : G4MagErrorStepper(EqRhs, numberOfVariables),
    fNumberOfVariables(numberOfVariables)
{
   dydxTemp = new G4double[fNumberOfVariables] ;
   yTemp    = new G4double[fNumberOfVariables] ;
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

  return ;
}  
