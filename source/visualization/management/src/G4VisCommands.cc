//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisCommands.cc,v 1.5 2001-08-17 23:02:17 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/ top level commands - John Allison  5th February 2001

#include "G4VisCommands.hh"

#include "G4VisManager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

////////////// /vis/enable ///////////////////////////////////////

G4VisCommandEnable::G4VisCommandEnable () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithABool("/vis/enable", this);
  fpCommand -> SetGuidance("/vis/enable [true|false]");
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

G4String G4VisCommandEnable::GetCurrentValue (G4UIcommand* command) {
  return G4String();
}

void G4VisCommandEnable::SetNewValue (G4UIcommand* command,
				      G4String newValue) {
  if (command == fpCommand) {
    G4bool enable (GetNewBoolValue(newValue));
    if (enable) fpVisManager->Enable();  // Printing is in vis manager.
    else fpVisManager->Disable();        // Printing is in vis manager.
  } else fpVisManager->Disable();        // Printing is in vis manager.
  // Note: Printing is in vis manager.
}

////////////// /vis/verbose ///////////////////////////////////////

G4VisCommandVerbose::G4VisCommandVerbose () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithAString("/vis/verbose", this);
  fpCommand -> SetGuidance(G4VisManager::VerbosityGuidanceString);
  fpCommand -> SetParameterName("verbosity", omitable=true);
  fpCommand -> SetDefaultValue("warnings");
}

G4VisCommandVerbose::~G4VisCommandVerbose () {
  delete fpCommand;
}

G4String G4VisCommandVerbose::GetCurrentValue (G4UIcommand* command) {
  return G4String();
}

void G4VisCommandVerbose::SetNewValue (G4UIcommand* command,
				       G4String newValue) {
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(newValue);
  fpVisManager->SetVerboseLevel(verbosity);
  // Always prints whatever the verbosity...
  G4cout << "Visualization verbosity changed to "
	 << G4VisManager::VerbosityString(verbosity) << G4endl;
}
