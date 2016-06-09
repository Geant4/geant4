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
// $Id: G4EvManMessenger.cc,v 1.5 2006/06/29 18:09:35 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

