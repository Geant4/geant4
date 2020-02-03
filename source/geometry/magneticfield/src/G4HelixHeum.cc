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
// G4HelixHeum implementation
//
// Simple Heum:
//        x_1 = x_0 + h *
//                1/4 * dx(t0,x0)  +
//                3/4 * dx(t0+2/3*h, x0+2/3*h*(dx(t0+h/3,x0+h/3*dx(t0,x0)))) 
//
//  Third order solver.
//
// Author: W.Wander <wwc@mit.edu>, 03/11/1998
// -------------------------------------------------------------------

#include "G4HelixHeum.hh"
#include "G4ThreeVector.hh"

G4HelixHeum::G4HelixHeum(G4Mag_EqRhs* EqRhs)
  : G4MagHelicalStepper(EqRhs)
{
}

G4HelixHeum::~G4HelixHeum()
{
}
  
void
G4HelixHeum::DumbStepper( const G4double  yIn[],
                          G4ThreeVector   Bfld,
                          G4double        h,
                          G4double        yOut[])
{
  const G4int nvar = 6 ;

  G4ThreeVector Bfield_Temp, Bfield_Temp2;
  G4double yTemp[6], yAdd1[6], yAdd2[6] , yTemp2[6];

  AdvanceHelix( yIn, Bfld, h, yAdd1 );
  
  AdvanceHelix( yIn, Bfld, h/3.0, yTemp );
  MagFieldEvaluate(yTemp,Bfield_Temp);

  AdvanceHelix( yIn, Bfield_Temp, (2.0 / 3.0) * h, yTemp2 );
  
  MagFieldEvaluate(yTemp2,Bfield_Temp2);

  AdvanceHelix( yIn, Bfield_Temp2, h, yAdd2 );

  for( G4int i = 0; i < nvar; ++i )
  {
    yOut[i] = ( 0.25 * yAdd1[i] + 0.75 * yAdd2[i]);
  }

  // NormaliseTangentVector( yOut );           
}  
