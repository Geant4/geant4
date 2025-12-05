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


#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::
PrimaryGeneratorMessenger(PrimaryGeneratorAction* gen)
  : fGen(gen)
{
  fBeamDir = new G4UIdirectory("/my_beam/");
  fBeamDir->SetGuidance("Example-specific commands to set beam properties");

  fKinECmd = new G4UIcmdWithADoubleAndUnit("/my_beam/kinE",this);
  fKinECmd->SetGuidance("Set the beam kinetic energy");
  fKinECmd->SetParameterName("kinE",false);
  fKinECmd->SetUnitCategory("Energy");
  fKinECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDECmd = new G4UIcmdWithADoubleAndUnit("/my_beam/DE",this);
  fDECmd->SetGuidance("Set the beam energy half-width, flat distribution");
  fDECmd->SetParameterName("DE",false);
  fDECmd->SetUnitCategory("Energy");
  fDECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fX0Cmd = new G4UIcmdWithADoubleAndUnit("/my_beam/X0",this);
  fX0Cmd->SetGuidance("Set X position of the center of the beam.");
  fX0Cmd->SetParameterName("X0",false);
  fX0Cmd->SetUnitCategory("Length");
  fX0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fY0Cmd = new G4UIcmdWithADoubleAndUnit("/my_beam/Y0",this);
  fY0Cmd->SetGuidance("Set Y position of the center of the beam.");
  fY0Cmd->SetParameterName("Y0",false);
  fY0Cmd->SetUnitCategory("Length");
  fY0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fZ0Cmd = new G4UIcmdWithADoubleAndUnit("/my_beam/Z0",this);
  fZ0Cmd->SetGuidance("Set Z position of the center of the beam.");
  fZ0Cmd->SetParameterName("Z0",false);
  fZ0Cmd->SetUnitCategory("Length");
  fZ0Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDXCmd = new G4UIcmdWithADoubleAndUnit("/my_beam/DX",this);
  fDXCmd->SetGuidance("Set the beam half-width for X, flat distribution");
  fDXCmd->SetParameterName("DX",false);
  fDXCmd->SetUnitCategory("Length");
  fDXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDYCmd = new G4UIcmdWithADoubleAndUnit("/my_beam/DY",this);
  fDYCmd->SetGuidance("Set the beam half-width for Y, flat distribution");
  fDYCmd->SetParameterName("DY",false);
  fDYCmd->SetUnitCategory("Length");
  fDYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fDZCmd = new G4UIcmdWithADoubleAndUnit("/my_beam/DZ",this);
  fDZCmd->SetGuidance("Set the beam half-width for Z, flat distribution");
  fDZCmd->SetParameterName("DZ",false);
  fDZCmd->SetUnitCategory("Length");
  fDZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fVerboseCmd = new G4UIcmdWithAnInteger("/my_beam/verbose", this);
  fVerboseCmd->SetGuidance("Set primary generator verbose");
  fVerboseCmd->SetParameterName("verb",false);
  fVerboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fBeamDir;
  delete fKinECmd;
  delete fDECmd;
  delete fX0Cmd;
  delete fY0Cmd;
  delete fZ0Cmd;
  delete fDXCmd;
  delete fDYCmd;
  delete fDZCmd;
  delete fVerboseCmd;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  if (command == fKinECmd)
    fGen->SetKinE(fKinECmd->GetNewDoubleValue(newValue));
  else if (command == fDECmd)
    fGen->SetDE(fDECmd->GetNewDoubleValue(newValue));
  else if (command == fX0Cmd)
    fGen->SetX0(fX0Cmd->GetNewDoubleValue(newValue));
  else if (command == fY0Cmd)
    fGen->SetY0(fY0Cmd->GetNewDoubleValue(newValue));
  else if (command == fZ0Cmd)
    fGen->SetZ0(fZ0Cmd->GetNewDoubleValue(newValue));
  else if (command == fDXCmd)
    fGen->SetDX(fDXCmd->GetNewDoubleValue(newValue));
  else if (command == fDYCmd)
    fGen->SetDY(fDYCmd->GetNewDoubleValue(newValue));
  else if (command == fDZCmd)
    fGen->SetDZ(fDZCmd->GetNewDoubleValue(newValue));
  else if (command == fVerboseCmd)
    fGen->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));
}

