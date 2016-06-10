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
// $Id: G4UCNBoundaryProcessMessenger.cc 69576 2013-05-08 13:48:13Z gcosmo $
//
///////////////////////////////////////////////////////////////////////
// UCN BoundaryProcess Messenger Class Implementation
///////////////////////////////////////////////////////////////////////
//
// File:        G4UCNBoundaryProcessMessenger.cc
// Description: Messenger class for UCN shutters
// Version:     1.0
// Created:     2014-05-12
// Author:      Peter Gumplinger
//              adopted from Geant4UCN by Peter Fierlinger (9.9.04) and
// Updated:
//
// mail:        gum@triumf.ca
//
///////////////////////////////////////////////////////////////////////

#include "G4UCNBoundaryProcess.hh"
#include "G4UCNBoundaryProcessMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UCNBoundaryProcessMessenger::
G4UCNBoundaryProcessMessenger(G4UCNBoundaryProcess* process)
  : theUCNBoundaryProcess (process)
{
  //G4cout << "Create messenger for the UCN boundary process " << G4endl;

  boundaryDir = new G4UIdirectory("/ucnboundary/");
  boundaryDir->SetGuidance("savetofile parameters");

  VerboseCmd = new G4UIcmdWithAnInteger("/ucnboundary/verbose",this);
  VerboseCmd->SetGuidance("Set verbose level" );
  VerboseCmd->SetParameterName("level",true);
  VerboseCmd->SetDefaultValue(1);
  VerboseCmd->AvailableForStates(G4State_Idle,G4State_PreInit);

  MicroRoughnessCmd =
                new G4UIcmdWithABool("/ucnboundary/MicroRoughness",this);
  MicroRoughnessCmd->
                     SetGuidance("Decide if MicroRoughness Model is activated");
  MicroRoughnessCmd->SetParameterName("MicroRough",false);
  MicroRoughnessCmd->SetDefaultValue(true);
  MicroRoughnessCmd->AvailableForStates(G4State_Idle,G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UCNBoundaryProcessMessenger::~G4UCNBoundaryProcessMessenger()
{
  if (VerboseCmd) delete VerboseCmd;
  if (MicroRoughnessCmd) delete MicroRoughnessCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UCNBoundaryProcessMessenger::SetNewValue(G4UIcommand* command,
                                                G4String newValue)
{
  if( command == VerboseCmd ) theUCNBoundaryProcess->
                        SetVerboseLevel(VerboseCmd->GetNewIntValue(newValue));
  if( command == MicroRoughnessCmd ) theUCNBoundaryProcess->
              SetMicroRoughness(MicroRoughnessCmd->GetNewBoolValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
