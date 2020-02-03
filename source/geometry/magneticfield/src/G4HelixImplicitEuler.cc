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
// G4HelixImplicitEuler implementation
//
//  Helix Implicit Euler:
//        x_1 = x_0 + 1/2 * ( helix(h,t_0,x_0)
//                          + helix(h,t_0+h,x_0+helix(h,t0,x0) ) )
//  Second order solver.
//  Take the current derivative and add it to the current position.
//  Take the output and its derivative. Add the mean of both derivatives
//  to form the final output
//
// Author: W.Wander <wwc@mit.edu>, 03/11/1998
// -------------------------------------------------------------------------

#include "G4HelixImplicitEuler.hh"
#include "G4ThreeVector.hh"

G4HelixImplicitEuler::G4HelixImplicitEuler(G4Mag_EqRhs *EqRhs)
  : G4MagHelicalStepper(EqRhs)
{
}

G4HelixImplicitEuler::~G4HelixImplicitEuler()
{
}
  
void
G4HelixImplicitEuler::DumbStepper( const G4double yIn[],
                                   G4ThreeVector  Bfld,
                                   G4double       h,
                                   G4double       yOut[])
{
  const G4int nvar = 6 ;
  G4double yTemp[6], yTemp2[6];
  G4ThreeVector Bfld_endpoint;

  // Step forward like in the explicit euler case
  //
  AdvanceHelix( yIn, Bfld, h, yTemp);

  // now obtain the new field value at the new point
  //
  MagFieldEvaluate(yTemp, Bfld_endpoint);      

  // and also advance along a helix for this field value
  //
  AdvanceHelix( yIn, Bfld_endpoint, h, yTemp2);

  // we take the average
  //
  for( G4int i = 0; i < nvar; ++i )
  {
    yOut[i] = 0.5 * ( yTemp[i] + yTemp2[i] );
  }

  // NormaliseTangentVector( yOut );           
}  
