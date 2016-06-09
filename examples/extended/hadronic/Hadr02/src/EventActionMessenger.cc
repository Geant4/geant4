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
/// \file hadronic/Hadr02/src/EventActionMessenger.cc
/// \brief Implementation of the EventActionMessenger class
//
// $Id$
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
#include "HistoManager.hh"

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

  ionCmd = new G4UIcmdWithAString("/testhadr/ionPhysics", this);
  ionCmd->SetGuidance("Added ion physics");
  ionCmd->SetGuidance("  Choice : FTF DPMJET");
  ionCmd->SetParameterName("ion",true);
  ionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
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
  delete ionCmd;
  delete printCmd;   
  delete dCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == drawCmd)
    {eventAction->SetDrawFlag(newValue);}

  if(command == ionCmd)
    {HistoManager::GetPointer()->SetIonPhysics(newValue);}
    
  if(command == printCmd)
    {eventAction->SetPrintModulo(printCmd->GetNewIntValue(newValue));}           

  if(command == dCmd)
    {eventAction->AddEventToDebug(dCmd->GetNewIntValue(newValue));}           
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
