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

G4String G4GAGTreeMessenger::GetCurrentValue(G4UIcommand*) {
  return "0";
}

void G4GAGTreeMessenger::SetNewValue(G4UIcommand*, G4String newValue) {
  fpGAGTree->SetVerbosity
    (fpCommandVerbose->GetNewIntValue(newValue));
  G4cout << "G4GAGTree verbosity now "
	 << fpGAGTree->GetVerbosity()
	 << G4endl;  
}
