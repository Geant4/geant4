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
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
//    *********************************
//    *                               *
//    *    RemSimRunMessenger.cc      *
//    *                               *
//    *********************************
//
//
// $Id: RemSimRunMessenger.cc,v 1.2 2004/05/22 12:57:07 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 

#include "RemSimRunMessenger.hh"
#include "RemSimRunAction.hh"
#include "RemSimRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

RemSimRunMessenger::RemSimRunMessenger( RemSimRunAction* arun): run(arun)
{ 
  runDir = new G4UIdirectory("/run/");
  runDir->SetGuidance("Select data.txt.");
        
  readFileCmd = new G4UIcmdWithAString("/run/data",this);
  readFileCmd->SetGuidance("Select data.txt/ write the name of the file"); 
  readFileCmd->SetParameterName("choice",true);
  readFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);         
 }

RemSimRunMessenger::~RemSimRunMessenger()
{
  delete readFileCmd;
  delete runDir;
}

void RemSimRunMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
   if(command == readFileCmd) run -> Read(newValue);
}

