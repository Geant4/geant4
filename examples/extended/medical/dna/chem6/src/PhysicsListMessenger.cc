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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
  : G4UImessenger()
  , fPhysicsList(pPhys)
{
  fPhysDir = std::make_unique<G4UIdirectory>("/chem6/",false);
  fPhysDir->SetGuidance("Time Step Model commands");
  fPhysDir->AvailableForStates(G4State_PreInit);

  fListCmd = std::make_unique<G4UIcmdWithAString>("/chem6/TimeStepModel", this);
  fListCmd->SetGuidance("Add modula chemistry list.");
  fListCmd->SetParameterName("ChemList", false);
  fListCmd->SetCandidates("SBS IRT IRT_syn");
  fListCmd->SetToBeBroadcasted(false);
  fListCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fListCmd.get())
  {
    if(newValue == "SBS")
    {
      fPhysicsList->RegisterTimeStepModel(TimeStepModel::fSBS);
    }
    else if(newValue == "IRT")
    {
      fPhysicsList->RegisterTimeStepModel(TimeStepModel::fIRT);
    }
    else if(newValue == "IRT_syn")
    {
      fPhysicsList->RegisterTimeStepModel(TimeStepModel::fIRT_syn);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
