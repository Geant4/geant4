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
// $Id: G4VisCommandTemplates.cc,v 1.4.4.1 2001/06/28 19:16:14 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
// Visualization Command Messenger Template private functions.
// John Allison  8th September 1997

#include "G4VisCommandTemplates.hh"

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4ios.hh"

// OLD STYLE!!
G4UIcommand* G4VisButtonCommandMessengerRegister (G4String commandName,
						  G4String guidance,
						  G4UImessenger* messenger) {
  G4UIcommand* command;
  G4UIparameter* param;
  command = new G4UIcommand (commandName, messenger);
  command -> SetGuidance (guidance);
  param   =  new G4UIparameter ("Selector", 'c', false);
  param   -> SetDefaultValue  ("?");
  command -> SetParameter     (param);
  return command;
}

// OLD STYLE!!
G4int G4VisButtonCommandMessengerInterpret (G4String newValues) {
  G4String& choice = newValues;
  G4int iSelector = -1;
  if (choice.compareTo ("0") == 0 ||
      choice.compareTo ("off", G4String::ignoreCase) == 0 ||
      choice.compareTo ("false", G4String::ignoreCase) == 0) iSelector = 0;
  if (choice.compareTo ("1") == 0 ||
      choice.compareTo ("on",  G4String::ignoreCase) == 0 ||
      choice.compareTo ("true", G4String::ignoreCase) == 0) iSelector = 1;
  if (iSelector < 0 || choice.isNull ()) {
    G4cout << "Choice not recognised."
      "\n  Choice is 0 (or off or false) or 1 (or on or true)." << G4endl;;
  }
  return iSelector;
}
