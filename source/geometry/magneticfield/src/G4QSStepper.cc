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
// G4QSStepper implementation
//
// Authors: Lucio Santi, Rodrigo Castro - 2018-2021.
// --------------------------------------------------------------------

#include "G4QSStepper.hh"

template<>
G4QSStepper_QSS3::G4QSStepper(G4EquationOfMotion *EqRhs,
                              G4int numberOfVariables,
                              G4bool primary)
  : G4QSStepper(new G4QSS3(G4QSStepper_QSS3::build_simulator()),
                EqRhs, numberOfVariables, primary)
{
}

template<>
G4QSStepper_QSS2::G4QSStepper(G4EquationOfMotion *EqRhs,
                              G4int numberOfVariables,
                              G4bool primary)
  : G4QSStepper(new G4QSS2(G4QSStepper_QSS2::build_simulator()),
                EqRhs, numberOfVariables, primary)
{
}

template<>
G4QSStepper_QSS2 *G4QSStepper_QSS2::build_QSS2(G4EquationOfMotion *EqRhs,
                                               G4int noIntegrationVariables,
                                               G4bool primary)
{
  return new G4QSStepper<G4QSS2>(EqRhs, noIntegrationVariables, primary);
}

template<>
G4QSStepper_QSS3 *G4QSStepper_QSS3::build_QSS3(G4EquationOfMotion *EqRhs,
                                               G4int noIntegrationVariables,
                                               G4bool primary)
{
  return new G4QSStepper<G4QSS3>(EqRhs, noIntegrationVariables, primary);
}
