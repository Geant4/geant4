// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITreeMessenger.cc,v 1.3 2001-06-15 07:12:37 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#include "G4ASCIITreeMessenger.hh"

#include "G4ASCIITree.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

G4ASCIITreeMessenger::G4ASCIITreeMessenger
(G4ASCIITree* ASCIITree):
  fpASCIITree(ASCIITree) {
  G4bool omitable, currentAsDefault;
  fpDirectory = new G4UIdirectory ("/vis/ASCIITree/");
  fpDirectory -> SetGuidance ("Commands for ASCIITree control.");
  fpCommandVerbose = new G4UIcmdWithAnInteger ("/vis/ASCIITree/verbose", this);
  fpCommandVerbose -> SetGuidance ("/vis/ASCIITree/verbose [<verbosity>]");
  fpCommandVerbose -> SetGuidance
    ("0 (default) mimimum - 10 maximum printing.");
  fpCommandVerbose -> SetParameterName ("verbosity",
					omitable = true,
					currentAsDefault = true);
}

G4ASCIITreeMessenger::~G4ASCIITreeMessenger() {
  delete fpCommandVerbose;
  delete fpDirectory;
}

G4String G4ASCIITreeMessenger::GetCurrentValue
(G4UIcommand* command) {
  return "0";
}

void G4ASCIITreeMessenger::SetNewValue
(G4UIcommand* command, G4String newValue) {
  fpASCIITree->SetVerbosity
    (fpCommandVerbose->GetNewIntValue(newValue));
  G4cout << "G4ASCIITree verbosity now "
	 << fpASCIITree->GetVerbosity()
	 << G4endl;  
}
