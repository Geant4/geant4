// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayXMessenger.cc,v 2.4 1998/07/16 01:57:32 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// GEANT4 Ray tracing for X windows Messenger - John Allison 5th June 1997

#ifdef G4VIS_BUILD_RAYX_DRIVER

#include "G4RayXMessenger.hh"

#include "G4RayX.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

G4RayXMessenger::G4RayXMessenger (G4RayX* pRayX):
fpRayX (pRayX)
{
  G4UIcommand* command;
  G4UIparameter* param;

  /////////////////////////////////////////////////  /vis~/ray/X...  ////
  //vis \hline
  //vis /vis~/ray/X... &&
  //vis ...any X-special commands for ray tracing... \\%
  //command = new G4UIcommand ("/vis~/ray/X...", this);
  //command -> SetGuidance ("...any X-special commands for ray tracing...");
}

G4RayXMessenger::~G4RayXMessenger () {}

void G4RayXMessenger::SetNewValue (G4UIcommand* command, G4String newValues)
{
  G4RayMessenger::SetNewValue (command, newValues);
}

G4String G4RayXMessenger::GetCurrentValue (G4UIcommand* command) {
  return "";
}

#endif
