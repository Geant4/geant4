// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EvManMessenger.cc,v 1.2 1999-12-15 14:49:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------

#include "G4EvManMessenger.hh"
#include "G4EventManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

G4EvManMessenger::G4EvManMessenger(G4EventManager * fEvMan)
:fEvManager(fEvMan)
{
  eventDirectory = new G4UIdirectory("/event/");
  eventDirectory->SetGuidance("EventManager control commands.");

  abortCmd = new G4UIcmdWithoutParameter("/event/abort",this);
  abortCmd->SetGuidance("Abort current event.");
  abortCmd->AvailableForStates(EventProc);

  verboseCmd = new G4UIcmdWithAnInteger("/event/verbose",this);
  verboseCmd->SetGuidance("Set Verbose level of event management category.");
  verboseCmd->SetGuidance(" 0 : Silent");
  verboseCmd->SetGuidance(" 1 : Stacking information");
  verboseCmd->SetGuidance(" 2 : More...");
  verboseCmd->SetParameterName("level",false);
  verboseCmd->SetRange("level>=0");
}

G4EvManMessenger::~G4EvManMessenger()
{
  delete abortCmd;
  delete verboseCmd;
  delete eventDirectory;
}

void G4EvManMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == verboseCmd )
  { fEvManager->SetVerboseLevel(verboseCmd->GetNewIntValue(newValues)); }
  if( command == abortCmd )
  { fEvManager->AbortCurrentEvent(); }
}

G4String G4EvManMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command == verboseCmd )
  { cv = verboseCmd->ConvertToString(fEvManager->GetVerboseLevel()); }
  return cv;
}

