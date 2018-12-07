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
/// \file RE05/src/RE05ActionInitialization.cc
/// \brief Implementation of the RE05ActionInitialization class
//

#include "RE05ActionInitialization.hh"
#include "RE05PrimaryGeneratorAction.hh"
#include "RE05RunAction.hh"
#include "RE05EventAction.hh"
#include "RE05StackingAction.hh"
#include "RE05TrackingAction.hh"
#include "RE05SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05ActionInitialization::RE05ActionInitialization()
: G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05ActionInitialization::~RE05ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05ActionInitialization::Build() const
{
  G4VUserPrimaryGeneratorAction* gen_action = new RE05PrimaryGeneratorAction;
  SetUserAction(gen_action);
  //
  G4UserRunAction* run_action = new RE05RunAction;
  SetUserAction(run_action);
  //
  G4UserEventAction* event_action = new RE05EventAction;
  SetUserAction(event_action);
  //
  G4UserStackingAction* stacking_action = new RE05StackingAction;
  SetUserAction(stacking_action);
  //
  G4UserTrackingAction* tracking_action = new RE05TrackingAction;
  SetUserAction(tracking_action);
  //
  G4UserSteppingAction* stepping_action = new RE05SteppingAction;
  SetUserAction(stepping_action);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05ActionInitialization::BuildForMaster() const
{
  //
  G4UserRunAction* run_action = new RE05RunAction;
  SetUserAction(run_action);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
