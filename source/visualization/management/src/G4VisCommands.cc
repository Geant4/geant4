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
// $Id: G4VisCommands.cc,v 1.14 2006/06/29 21:29:34 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $

// /vis/ top level commands - John Allison  5th February 2001

#include "G4VisCommands.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

////////////// /vis/enable ///////////////////////////////////////

G4VisCommandEnable::G4VisCommandEnable () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithABool("/vis/enable", this);
  fpCommand -> SetGuidance("Enables/disables visualization system.");
  fpCommand -> SetParameterName("enabled", omitable=true);
  fpCommand -> SetDefaultValue(true);

  fpCommand1 = new G4UIcmdWithoutParameter("/vis/disable", this);
  fpCommand1 -> SetGuidance("Disables visualization system.");
}

G4VisCommandEnable::~G4VisCommandEnable () {
  delete fpCommand;
  delete fpCommand1;
}

G4String G4VisCommandEnable::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandEnable::SetNewValue (G4UIcommand* command,
				      G4String newValue) {
  if (command == fpCommand) {
    G4bool enable = G4UIcommand::ConvertToBool(newValue);
    if (enable) fpVisManager->Enable();  // Printing is in vis manager.
    else fpVisManager->Disable();        // Printing is in vis manager.
  } else fpVisManager->Disable();        // Printing is in vis manager.
  // Note: Printing is in vis manager.
}

////////////// /vis/list ///////////////////////////////////////

G4VisCommandList::G4VisCommandList ()
{
  G4bool omitable;

  fpCommand = new G4UIcmdWithAString("/vis/list", this);
  fpCommand -> SetGuidance("Lists visualization parameters.");
  fpCommand -> SetParameterName("verbosity", omitable=true);
  fpCommand -> SetDefaultValue("warnings");
}

G4VisCommandList::~G4VisCommandList ()
{
  delete fpCommand;
}

G4String G4VisCommandList::GetCurrentValue (G4UIcommand*)
{
  return "";
}

void G4VisCommandList::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4String& verbosityString = newValue;
  //G4VisManager::Verbosity verbosity =
  //  fpVisManager->GetVerbosityValue(verbosityString);

  fpVisManager->PrintAvailableGraphicsSystems();
  G4cout << G4endl;
  fpVisManager->PrintAvailableModels();
  G4cout << G4endl;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand(G4String("/vis/viewer/list ! ") + verbosityString);
}

////////////// /vis/verbose ///////////////////////////////////////

G4VisCommandVerbose::G4VisCommandVerbose () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithAString("/vis/verbose", this);
  for (size_t i = 0; i < G4VisManager::VerbosityGuidanceStrings.size(); ++i) {
    fpCommand -> SetGuidance(G4VisManager::VerbosityGuidanceStrings[i]);
  }
  fpCommand -> SetParameterName("verbosity", omitable=true);
  fpCommand -> SetDefaultValue("warnings");
}

G4VisCommandVerbose::~G4VisCommandVerbose () {
  delete fpCommand;
}

G4String G4VisCommandVerbose::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandVerbose::SetNewValue (G4UIcommand*,
				       G4String newValue) {
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(newValue);
  fpVisManager->SetVerboseLevel(verbosity);
  // Always prints whatever the verbosity...
  G4cout << "Visualization verbosity changed to "
	 << G4VisManager::VerbosityString(verbosity) << G4endl;
}
