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
// $Id: TstVAEventActionMessenger.cc,v 1.5 2002-12-04 19:09:31 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#include "globals.hh"

#include "TstVAEventActionMessenger.hh"

#include "TstVAEventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

TstVAEventActionMessenger::TstVAEventActionMessenger(TstVAEventAction* EvAct)
:eventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none, charged(default),neutral, all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged neutral all");
  DrawCmd->AvailableForStates(G4State_Idle);
  
  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);     
}

TstVAEventActionMessenger::~TstVAEventActionMessenger()
{
  delete DrawCmd;
  delete PrintCmd;   
}

void TstVAEventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
       
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}              
}
