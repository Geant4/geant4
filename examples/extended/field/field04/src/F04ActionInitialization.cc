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
/// \file F04ActionInitialization.cc
/// \brief Implementation of the F04ActionInitialization class

#include "F04ActionInitialization.hh"
#include "F04PrimaryGeneratorAction.hh"
#include "F04RunAction.hh"
#include "F04EventAction.hh"
#include "F04TrackingAction.hh"
#include "F04StackingAction.hh"
#include "F04SteppingAction.hh"
#include "F04SteppingVerbose.hh"

#include "F04DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04ActionInitialization::F04ActionInitialization
                            (F04DetectorConstruction* detConstruction)
 : G4VUserActionInitialization(),
   fDetConstruction(detConstruction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04ActionInitialization::~F04ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04ActionInitialization::BuildForMaster() const
{
  SetUserAction(new F04RunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04ActionInitialization::Build() const
{
  SetUserAction(new F04PrimaryGeneratorAction(fDetConstruction));

  F04RunAction* runAction = new F04RunAction();
  SetUserAction(runAction);
  F04EventAction* eventAction = new F04EventAction(runAction);
  SetUserAction(eventAction);
  SetUserAction(new F04TrackingAction());
  SetUserAction(new F04StackingAction());
  SetUserAction(new F04SteppingAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose* F04ActionInitialization::InitializeSteppingVerbose() const
{
  return new F04SteppingVerbose();
}
