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
/// \file optical/OpNovice2/src/SteppingMessenger.cc
/// \brief Implementation of the SteppingMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingMessenger.hh"
#include "SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingMessenger::SteppingMessenger(SteppingAction* steppingAction)
  : G4UImessenger(),
    fSteppingAction(steppingAction)
{
  fSteppingDir = new G4UIdirectory("/opnovice2/stepping/");
  fSteppingDir->SetGuidance("Stepping control");

  fKillOnSecondSurfaceCmd =
    new G4UIcmdWithABool("/opnovice2/stepping/killOnSecondSurface", this);
  fKillOnSecondSurfaceCmd->SetGuidance(
    "Kill the optical photon when it reaches a second surface. "
    "Useful for visualizing boundary scattering.");
  fKillOnSecondSurfaceCmd->SetDefaultValue(false);
  fKillOnSecondSurfaceCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingMessenger::~SteppingMessenger()
{
  delete fSteppingDir;
  delete fKillOnSecondSurfaceCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingMessenger::SetNewValue(G4UIcommand* command,
                                    G4String newValue)
{
  if(command == fKillOnSecondSurfaceCmd)
  {
    fSteppingAction->SetKillOnSecondSurface(
      G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
