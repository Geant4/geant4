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
// $Id: EventActionMessenger.cc,v 1.1 2003-05-27 13:44:49 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "EventActionMessenger.hh"

#include "EventAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"


EventActionMessenger::EventActionMessenger(EventAction* EvAct)
 :eventAction(EvAct)
{
  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);     
}


EventActionMessenger::~EventActionMessenger()
{
  delete PrintCmd;   
}


void EventActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}
}


