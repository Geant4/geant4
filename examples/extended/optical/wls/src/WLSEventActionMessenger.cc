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
/// \file optical/wls/src/WLSEventActionMessenger.cc
/// \brief Implementation of the WLSEventActionMessenger class
//
//
//

#include "globals.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "WLSEventAction.hh"
#include "WLSEventActionMessenger.hh"

WLSEventActionMessenger::WLSEventActionMessenger(WLSEventAction* EvAct)
     :eventAction(EvAct)
{
  setVerboseCmd = new G4UIcmdWithAnInteger("/event/setverbose",this);
  setVerboseCmd->SetGuidance("Set verbose level ." );
  setVerboseCmd->SetParameterName("level",true);
  setVerboseCmd->SetDefaultValue(0);

  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none,charged, all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged all");
  DrawCmd->AvailableForStates(G4State_Idle);

  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);
}

WLSEventActionMessenger::~WLSEventActionMessenger()
{
  delete setVerboseCmd;
  delete DrawCmd;
  delete PrintCmd;
}

void WLSEventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if (command == setVerboseCmd)
    eventAction->SetEventVerbose(setVerboseCmd->GetNewIntValue(newValue));
 
  if (command == DrawCmd)
    eventAction->SetDrawFlag(newValue);

  if (command == PrintCmd)
    eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));
}
