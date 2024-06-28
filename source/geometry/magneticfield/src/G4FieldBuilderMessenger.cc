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

/// \file G4FieldBuilderMessenger.cc
/// \brief Implementation of the G4FieldBuilderMessenger class
///
/// \author I. Hrivnacova; IJCLab, Orsay

#include "G4FieldBuilderMessenger.hh"
#include "G4FieldBuilder.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

//_____________________________________________________________________________
G4FieldBuilderMessenger::G4FieldBuilderMessenger(G4FieldBuilder* fieldBuilder)
  : fFieldBuilder(fieldBuilder)
{
  // Standard constructor

  G4String directoryName = "/field/";
  fDirectory = new G4UIdirectory(directoryName);
  fDirectory->SetGuidance("Magnetic (or other type) field control commands.");

  G4String commandName = directoryName;
  commandName.append("verboseLevel");
  fVerboseLevelCmd = new G4UIcmdWithAnInteger(commandName, this);
  fVerboseLevelCmd->SetGuidance("Set verbose level");
  fVerboseLevelCmd->SetParameterName("VerboseLevel", false);
  // add possible values
  fVerboseLevelCmd->AvailableForStates(
    G4State_PreInit, G4State_Init, G4State_Idle);
}

//_____________________________________________________________________________
G4FieldBuilderMessenger::~G4FieldBuilderMessenger()
{
  // Destructor

  delete fDirectory;
  delete fVerboseLevelCmd;
}

//
// public methods
//

//_____________________________________________________________________________
void G4FieldBuilderMessenger::SetNewValue(
  G4UIcommand* command, G4String newValues)
{
  // Apply command to the associated object.

  if (command == fVerboseLevelCmd) {
    fFieldBuilder->SetVerboseLevel(fVerboseLevelCmd->GetNewIntValue(newValues));
    return;
  }
}
