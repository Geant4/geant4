// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4WoMessenger.cc,v 1.2.8.1 1999/12/07 20:49:06 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#ifdef G4UI_BUILD_WO_SESSION

#include "G4WoMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

#include <Wo.h>

G4WoMessenger::G4WoMessenger(OShell a_shell)
:shell(NULL)
,WoDirectory(NULL)
,sendExitCmd(NULL)
,oshCmd(NULL)
{
  shell                   = a_shell;

  WoDirectory = new G4UIdirectory("/Wo/");
  WoDirectory->SetGuidance  ("Wo control commands.");

  sendExitCmd = new G4UIcmdWithoutParameter("/Wo/sendExit",this);
  sendExitCmd->SetGuidance  ("Exit Wo event event loop..");

  oshCmd = new G4UIcmdWithAString("/Wo/osh",this);
  oshCmd->SetGuidance ("osh interpreter.");
  oshCmd->SetParameterName("Command",true);
  oshCmd->SetDefaultValue ("");
}
G4WoMessenger::~G4WoMessenger() {
  if(oshCmd!=NULL)      delete oshCmd;
  if(sendExitCmd!=NULL) delete sendExitCmd;
  if(WoDirectory!=NULL) delete WoDirectory;
}

void G4WoMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if(command==oshCmd)  {
    OShellExecute (shell,(char*)newValues.data());
  } else if(command==sendExitCmd)  {
    WoSendExit ();
  }
}
G4String G4WoMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  return cv;
}

#endif
