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
/// \file radiobiology/src/PhysicsListMessenger.cc
/// \brief Implementation of the RadioBio::PhysicsListMessenger class

#include "PhysicsListMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "PhysicsList.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* phys) : fPhysicsList(phys)
{
  // Directory for Physics commands
  fPhysDir = new G4UIdirectory("/Physics/");
  fPhysDir->SetGuidance("Commands to activate physics models and set cuts");

  // Add a Physics List
  fListCmd = new G4UIcmdWithAString("/Physics/addPhysics", this);
  fListCmd->SetGuidance("Add physics list.");
  fListCmd->SetGuidance("The available choices are: ");
  fListCmd->SetGuidance("standard_opt4: Only E.M. Physics");
  fListCmd->SetGuidance("HADRONTHERAPY_1: E.M. Physics (standard_opt4) + Hadron Physics HP");
  fListCmd->SetGuidance("HADRONTHERAPY_2: E.M. Physics (standard_opt4) + Hadron Physics");
  fListCmd->SetCandidates("standard_opt4 HADRONTHERAPY_1 HADRONTHERAPY_2");
  fListCmd->SetDefaultValue("HADRONTHERAPY_1");
  fListCmd->SetParameterName("PList", false);
  fListCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fPhysDir;
  delete fListCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fListCmd) {
    fPhysicsList->AddPhysicsList(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio