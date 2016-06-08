// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCameraMessengers.cc,v 1.2.8.1 1999/12/07 20:54:00 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// Messengers for /vis~/camera commands - John Allison  6th April 1998.

#include "G4VisCommandsCameraMessengers.hh"

#include "G4UIcommand.hh"

G4VisCommandsCameraSetMessenger::G4VisCommandsCameraSetMessenger () {
  fpCommand = new G4UIcommand ("/new_vis/camera/set", this);
  fpCommand -> SetGuidance
    ("Sets and updates view parameters of current view.");
}

G4VisCommandsCameraSetMessenger::~G4VisCommandsCameraSetMessenger () {
  delete fpCommand;
}

G4String G4VisCommandsCameraSetMessenger::GetCurrentValue
(G4UIcommand * command) {
  return "";
}

void G4VisCommandsCameraSetMessenger::SetNewValue (G4UIcommand* command, G4String newValues) {
}
