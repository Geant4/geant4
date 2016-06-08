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
// Satoshi Tanaka 31th May 2001
//
// A messenger for G4GAGTree driver.

#include "G4GAGTreeMessenger.hh"

#include "G4GAGTree.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

G4GAGTreeMessenger::G4GAGTreeMessenger
(G4GAGTree* GAGTree):
  fpGAGTree(GAGTree) {
  G4bool omitable, currentAsDefault;
  fpDirectory = new G4UIdirectory ("/vis/GAGTree/");
  fpDirectory -> SetGuidance ("Commands for GAGTree control.");
  fpCommandVerbose = new G4UIcmdWithAnInteger ("/vis/GAGTree/verbose", this);
  fpCommandVerbose -> SetGuidance ("/vis/GAGTree/verbose [<verbosity>]");
  fpCommandVerbose -> SetGuidance
    ("0 (default) mimimum - 10 maximum printing.");
  fpCommandVerbose -> SetParameterName ("verbosity",
					omitable = true,
					currentAsDefault = true);
}

G4GAGTreeMessenger::~G4GAGTreeMessenger() {
  delete fpCommandVerbose;
  delete fpDirectory;
}

G4String G4GAGTreeMessenger::GetCurrentValue
(G4UIcommand* command) {
  return "0";
}

void G4GAGTreeMessenger::SetNewValue
(G4UIcommand* command, G4String newValue) {
  fpGAGTree->SetVerbosity
    (fpCommandVerbose->GetNewIntValue(newValue));
  G4cout << "G4GAGTree verbosity now "
	 << fpGAGTree->GetVerbosity()
	 << G4endl;  
}
