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
// $Id: EventActionMessenger.cc,v 1.2 2006-06-06 19:48:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// EventActionMessenger
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

#include "EventActionMessenger.hh"

#include "EventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventActionMessenger::EventActionMessenger(EventAction* EvAct)
:eventAction(EvAct)
{ 
  drawCmd = new G4UIcmdWithAString("/testhadr/DrawTracks", this);
  drawCmd->SetGuidance("Draw the tracks in the event");
  drawCmd->SetGuidance("  Choice : neutral, charged, all");
  drawCmd->SetParameterName("choice",true);
  drawCmd->SetDefaultValue("all");
  drawCmd->SetCandidates("none charged all");
  drawCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  printCmd = new G4UIcmdWithAnInteger("/testhadr/PrintModulo",this);
  printCmd->SetGuidance("Print events modulo n");
  printCmd->SetParameterName("EventNb",false);
  printCmd->SetRange("EventNb>0");
  printCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      

  dCmd = new G4UIcmdWithAnInteger("/testhadr/DebugEvent",this);
  dCmd->SetGuidance("D event to debug");
  dCmd->SetParameterName("fNb",false);
  dCmd->SetRange("fNb>0");
  dCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventActionMessenger::~EventActionMessenger()
{
  delete drawCmd;
  delete printCmd;   
  delete dCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == drawCmd)
    {eventAction->SetDrawFlag(newValue);}
    
  if(command == printCmd)
    {eventAction->SetPrintModulo(printCmd->GetNewIntValue(newValue));}           

  if(command == dCmd)
    {eventAction->AddEventToDebug(dCmd->GetNewIntValue(newValue));}           
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
