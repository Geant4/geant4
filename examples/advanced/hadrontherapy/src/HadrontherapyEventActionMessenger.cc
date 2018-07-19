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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#include "HadrontherapyEventActionMessenger.hh"

#include "HadrontherapyEventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"


/////////////////////////////////////////////////////////////////////////////
HadrontherapyEventActionMessenger::HadrontherapyEventActionMessenger(HadrontherapyEventAction* EvAct)
:eventAction(EvAct)
{
  eventDir = new G4UIdirectory("/event/");
  eventDir->SetGuidance("Permits controls on simulation events");
 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none,charged, all, neutral");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged all neutral");
  DrawCmd->AvailableForStates(G4State_Idle);
  
  PrintCmd = new G4UIcmdWithAnInteger("/event/printEventNumber",this);
  PrintCmd->SetGuidance("Print the event number of modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);      
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyEventActionMessenger::~HadrontherapyEventActionMessenger()
{
  delete DrawCmd;
  delete PrintCmd;
  delete eventDir;         
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyEventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
    
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}           
 }
