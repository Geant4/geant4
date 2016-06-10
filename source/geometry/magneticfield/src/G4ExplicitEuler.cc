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
// $Id: G4ExplicitEuler.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
//
//  Explicit Euler: x_1 = x_0 + h * dx_0
//
//  most simple approach for solving linear differential equations.
//  Take the current derivative and add it to the current position.
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// -------------------------------------------------------------------

#include "G4ExplicitEuler.hh"
#include "G4ThreeVector.hh"

//////////////////////////////////////////////////////////////////////////
//
// Constructor

G4ExplicitEuler::G4ExplicitEuler(G4EquationOfMotion* EqRhs, 
                                 G4int numberOfVariables)
 : G4MagErrorStepper(EqRhs, numberOfVariables)
{
}


///////////////////////////////////////////////////////////////////////
//
// Destructor

G4ExplicitEuler::~G4ExplicitEuler()
{
}


///////////////////////////////////////////////////////////////////////
//
//

void
G4ExplicitEuler::DumbStepper( const G4double  yIn[],
			      const G4double  dydx[],
			            G4double  h,
			 	    G4double  yOut[]        )
{
  const G4int numberOfVariables= GetNumberOfVariables();

  // Initialise time to t0, needed when it is not updated by the integration.
  // yOut[7] = yIn[7];   //  Better to set it to NaN;  // TODO

  G4int i;

  for(i=0;i< numberOfVariables;i++)
  {
    yOut[i] = yIn[i] + h*dydx[i] ;             // 1st and only Step 
  }
  
  return ;
}  
