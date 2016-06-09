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
/// \file field/field06/src/F06Field.cc
/// \brief Implementation of the F06Field class
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F06Field.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4EqGravityField.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"

#include "G4ClassicalRK4.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06Field::F06Field() : G4UniformGravityField()
{
  fEquation = new G4EqGravityField(this);

  G4FieldManager* fFieldManager
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  fFieldManager->SetDetectorField(this);

  fStepper = new G4ClassicalRK4(fEquation,8);

  G4double minStep           = 0.01*mm;

  fChordFinder = new G4ChordFinder((G4MagneticField*)this,minStep,fStepper);
 
  // Set accuracy parameters
  G4double deltaChord        = 3.0*mm;
  fChordFinder->SetDeltaChord( deltaChord );

  G4double deltaOneStep      = 0.01*mm; 
  fFieldManager->SetAccuraciesWithDeltaOneStep(deltaOneStep);

  G4double deltaIntersection = 0.1*mm; 
  fFieldManager->SetDeltaIntersection(deltaIntersection);

  G4TransportationManager* fTransportManager =
         G4TransportationManager::GetTransportationManager();

  fieldPropagator = fTransportManager->GetPropagatorInField();

  G4double epsMin            = 2.5e-7*mm;
  G4double epsMax            = 0.05*mm;
 
  fieldPropagator->SetMinimumEpsilonStep(epsMin);
  fieldPropagator->SetMaximumEpsilonStep(epsMax);
 
  fFieldManager->SetChordFinder(fChordFinder);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F06Field::~F06Field()
{
  if (fEquation)    delete fEquation;
  if (fStepper)     delete fStepper;
  if (fChordFinder) delete fChordFinder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
