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
// $Id: G4EvManMessenger.cc,v 1.4 2002-12-04 18:43:15 asaim Exp $
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
  abortCmd->AvailableForStates(G4State_EventProc);

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

