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
// $Id: G4SimpleHeum.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
//  Simple Heum:
//        x_1 = x_0 + h *
//                1/4 * dx(t0,x0)  +
//                3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
//
// Third order solver.
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// -------------------------------------------------------------------

#include "G4SimpleHeum.hh"
#include "G4ThreeVector.hh"

///////////////////////////////////////////////////////////////////////////
//
// Constructor

G4SimpleHeum::G4SimpleHeum(G4EquationOfMotion *EqRhs, G4int num_variables): 
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
                                 G4double  h,
                                 G4double  yOut[])
{
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
    yOut[i] = yIn[i] + h * (0.25 * dydx[i] + 0.75 * dydxTemp2[i]);
  }
      
  if ( fNumberOfVariables == 12 ) { NormalisePolarizationVector( yOut ); }
}  
