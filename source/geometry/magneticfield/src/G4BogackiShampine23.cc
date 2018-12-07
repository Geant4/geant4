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
//  Bogacki-Shampine - 4 - 3(2) non-FSAL implementation by Somnath Banerjee
//  Supervision / code review: John Apostolakis
//
//  Somnath's work was sponsored by Google as part of the Google Summer of
//  Code 2015, as part of the CERN / SFT organisation.p
// ===================================================================
//
//  Implementation of the method proposed in the publication
//   “A 3(2) pair of Runge - Kutta formulas,”
//  by P. Bogacki and L. F. Shampine,
//  Appl. Math. Lett., vol. 2, no. 4, pp. 321–325, Jan. 1989.
// 
// First version: 20 May 2015
//
//  History
// -----------------------------
//  Created by Somnath Banerjee on 20 May 2015
///////////////////////////////////////////////////////////////////////////////


/*

This contains the stepper function of the G4BogackiShampine23 class

The Bogacki shampine method has the following Butcher's tableau

0  |
1/2|1/2
3/4|0	3/4
1  |2/9	1/3	4/9
-------------------
   |2/9	1/3	4/9	0
   |7/24 1/4 1/3 1/8

*/

#include "G4BogackiShampine23.hh"
#include "G4LineSection.hh"
#include "G4FieldUtils.hh"

using namespace field_utils;

G4BogackiShampine23::G4BogackiShampine23(G4EquationOfMotion* EqRhs,
				                                 G4int integrationVariables): 
  G4MagIntegratorStepper(EqRhs, integrationVariables)
{

  SetIntegrationOrder(3);
  SetFSAL(true);
}

void G4BogackiShampine23::makeStep(const G4double yInput[],
                                   const G4double dydx[],
                                   const G4double hstep,
                                   G4double yOutput[],
                                   G4double* dydxOutput,
                                   G4double* yError) const
{

  G4double yTemp[G4FieldTrack::ncompSVEC];
  for(G4int i = GetNumberOfVariables(); i < GetNumberOfStateVariables(); ++i) 
    yOutput[i] = yTemp[i] = yInput[i];
  

  G4double ak2[G4FieldTrack::ncompSVEC],
           ak3[G4FieldTrack::ncompSVEC];

 const G4double b21 = 0.5 ,
                b31 = 0., b32 = 3.0 / 4.0,
                b41 = 2.0 / 9.0, b42 = 1.0 / 3.0, b43 = 4.0 / 9.0;

 const G4double dc1 = b41 - 7.0 / 24.0,  dc2 = b42 - 1.0 / 4.0,
  				      dc3 = b43 - 1.0 / 3.0, dc4 = - 1.0 / 8.0;
 
    // RightHandSide(yInput, dydx);
    for(G4int i = 0; i < GetNumberOfVariables(); ++i)
      yTemp[i] = yInput[i] + b21 * hstep * dydx[i];
    
    RightHandSide(yTemp, ak2);
    for(G4int i = 0; i < GetNumberOfVariables(); ++i)
      yTemp[i] = yInput[i] + hstep * (b31 * dydx[i] + b32 * ak2[i]);

    RightHandSide(yTemp, ak3);
    for(G4int i = 0; i < GetNumberOfVariables(); ++i)
      yOutput[i] = yInput[i] + hstep * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    
    if (dydxOutput && yError) {
      RightHandSide(yOutput, dydxOutput);
      for(G4int i = 0; i < GetNumberOfVariables(); ++i)
        yError[i] = hstep * (dc1 * dydx[i] + dc2 * ak2[i] + 
                             dc3 * ak3[i] + dc4 * dydxOutput[i]);
      
  }
}

void G4BogackiShampine23::Stepper(const G4double yInput[],
                                  const G4double dydx[],
                                  G4double hstep,
                                  G4double yOutput[],
                                  G4double yError[])
{
  copy(fyIn, yInput);
  copy(fdydx, dydx);
  fhstep = hstep;

  makeStep(fyIn, fdydx, fhstep, fyOut, fdydxOut, yError);

  copy(yOutput, fyOut);
}

void G4BogackiShampine23::Stepper(const G4double yInput[],
                                  const G4double dydx[],
                                  G4double hstep,
                                  G4double yOutput[],
                                  G4double yError[],
                                  G4double dydxOutput[])
{
  copy(fyIn, yInput);
  copy(fdydx, dydx);
  fhstep = hstep;

  makeStep(fyIn, fdydx, fhstep, fyOut, fdydxOut, yError);

  copy(yOutput, fyOut);
  copy(dydxOutput, fdydxOut);
}

G4double G4BogackiShampine23::DistChord() const
{
  G4double yMid[G4FieldTrack::ncompSVEC];
  makeStep(fyIn, fdydx, fhstep / 2., yMid);

  const G4ThreeVector begin = makeVector(fyIn, Value3D::Position);
  const G4ThreeVector mid = makeVector(yMid, Value3D::Position);
  const G4ThreeVector end = makeVector(fyOut, Value3D::Position);

  return G4LineSection::Distline(mid, begin, end);
}
