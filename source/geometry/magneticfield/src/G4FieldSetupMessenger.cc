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
//------------------------------------------------
// The Geant4 Virtual Monte Carlo package
// Copyright (C) 2007 - 2014 Ivana Hrivnacova
// All rights reserved.
//
// For the licensing terms see geant4_vmc/LICENSE.
// Contact: root-vmc@cern.ch
//-------------------------------------------------

/// \file G4FieldSetupMessenger.cc
/// \brief Implementation of the G4FieldSetupMessenger class
///
/// \author I. Hrivnacova; IJCLab, Orsay

#include "G4FieldSetupMessenger.hh"
#include "G4FieldSetup.hh"
#include "G4LogicalVolume.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

//_____________________________________________________________________________
G4FieldSetupMessenger::G4FieldSetupMessenger(G4FieldSetup* fieldSetup)
  : fFieldSetup(fieldSetup)
{
  // Standard constructor

  G4String directoryName = "/field/";
  if (fFieldSetup->GetLogicalVolume() != nullptr) {
    directoryName.append(fFieldSetup->GetLogicalVolume()->GetName());
    directoryName.append("/");
  }

  G4String commandName = std::move(directoryName);
  commandName.append("update");
  fUpdateCmd = new G4UIcmdWithoutParameter(commandName, this);
  fUpdateCmd->SetGuidance("Update field setup.");
  fUpdateCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);
}

//_____________________________________________________________________________
G4FieldSetupMessenger::~G4FieldSetupMessenger()
{
  // Destructor

  delete fUpdateCmd;
}

//
// public methods
//

//_____________________________________________________________________________
void G4FieldSetupMessenger::SetNewValue(
  G4UIcommand* command, G4String /*newValues*/)
{
  // Apply command to the associated object.

  if (command == fUpdateCmd) {
    G4cout << "Execute update command" << G4endl;
    fFieldSetup->Update();
    return;
  }
}
