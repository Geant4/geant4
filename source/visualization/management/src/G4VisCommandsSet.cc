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
// $Id: G4VisCommands.cc,v 1.24 2009-03-09 12:42:00 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/set - John Allison  21st March 2012
// Set quantities for use in appropriate commands.

#include "G4VisCommandsSet.hh"

#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include <cctype>

////////////// /vis/set/textColour ////////////////////////////////////

G4VisCommandSetTextColour::G4VisCommandSetTextColour ()
{
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/set/textColour", this);
  fpCommand -> SetGuidance
    ("Defines colour and opacity for future \"/vis/scene/add/text\" commands.");
  fpCommand -> SetGuidance("Default: blue and opaque.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter ("red", 's', omitable = true);
  parameter->SetGuidance
    ("Red component or a string, e.g., \"blue\" (green and blue parameters are ignored).");
  parameter->SetDefaultValue ("0.");
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("green", 'd', omitable = true);
  parameter->SetDefaultValue (0.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("blue", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  fpCommand->SetParameter (parameter);
  parameter = new G4UIparameter ("alpha", 'd', omitable = true);
  parameter->SetDefaultValue (1.);
  parameter->SetGuidance ("Opacity");
  fpCommand->SetParameter (parameter);
}

G4VisCommandSetTextColour::~G4VisCommandSetTextColour ()
{
  delete fpCommand;
}

G4String G4VisCommandSetTextColour::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetTextColour::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4String redOrString;
  G4double green, blue, opacity;
  std::istringstream iss(newValue);
  iss >> redOrString >> green >> blue >> opacity;

  G4Colour colour(0,0,1,1);  // Default blue and opaque.
  if (std::isalpha(redOrString(0))) {
    if (!G4Colour::GetColour(redOrString, colour)) {
      if (verbosity >= G4VisManager::warnings) {
	G4cout << "WARNING: Text colour \"" << redOrString
	       << "\" not found.  Defaulting to blue and opaque."
	       << G4endl;
      }
    }
  } else {
    colour = G4Colour
      (G4UIcommand::ConvertToDouble(redOrString), green, blue);
  }
  // Add opacity
  fCurrentTextColour = G4Colour
    (colour.GetRed(), colour.GetGreen(), colour.GetBlue(), opacity);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Text colour (for future \"text\" commands) has been set to "
	   << fCurrentTextColour << '.'
	   << G4endl;
  }
}

////////////// /vis/set/textLayout ////////////////////////////////////

G4VisCommandSetTextLayout::G4VisCommandSetTextLayout ()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/set/textLayout", this);
  fpCommand -> SetGuidance
    ("Defines layout future \"/vis/scene/add/text\" commands.");
  fpCommand -> SetGuidance
    ("\"left\" for left justification to provided coordinate.");
  fpCommand -> SetGuidance
    ("\"centre\" or \"center\" for text centered on provided coordinate.");
  fpCommand -> SetGuidance
    ("\"right\" for right justification to provided coordinate.");
  fpCommand -> SetGuidance("Default: left.");
  fpCommand -> SetParameterName("layout", omitable = false);
  fpCommand -> SetCandidates ("left centre center right");
  fpCommand -> SetDefaultValue ("left");
}

G4VisCommandSetTextLayout::~G4VisCommandSetTextLayout ()
{
  delete fpCommand;
}

G4String G4VisCommandSetTextLayout::GetCurrentValue (G4UIcommand*)
{
  return G4String();
}

void G4VisCommandSetTextLayout::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4Text::Layout layout = G4Text::left;
  if (newValue == "left") layout = G4Text::left;
  else if (newValue == "centre" || newValue == "center")
    layout = G4Text::centre;
  else if (newValue == "right") layout = G4Text::right;

  fCurrentTextLayout = layout;

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  if (verbosity >= G4VisManager::confirmations) {
    G4cout << "Text layout (for future \"text\" commands) has been set to \""
	   << fCurrentTextLayout << "\"."
	   << G4endl;
  }
}
