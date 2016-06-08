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
// $Id: G4SDmessenger.cc,v 1.3 2001/07/11 09:58:44 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// ---------------------------------------------------------------------

#include "G4SDmessenger.hh"
#include "G4SDManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

G4SDmessenger::G4SDmessenger(G4SDManager* SDManager):fSDMan(SDManager)
{
  hitsDir = new G4UIdirectory("/hits/");
  hitsDir->SetGuidance("Sensitive detectors and Hits");

  listCmd = new G4UIcmdWithoutParameter("/hits/list",this);
  listCmd->SetGuidance("List sensitive detector tree.");

  activeCmd = new G4UIcmdWithAString("/hits/activate",this);
  activeCmd->SetGuidance("Activate sensitive detector(s).");
  activeCmd->SetParameterName("detector",true);
  activeCmd->SetDefaultValue("/");

  inactiveCmd = new G4UIcmdWithAString("/hits/inactivate",this);
  inactiveCmd->SetGuidance("Inactivate sensitive detector(s).");
  inactiveCmd->SetParameterName("detector",true);
  inactiveCmd->SetDefaultValue("/");

  verboseCmd = new G4UIcmdWithAnInteger("/hits/verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level.");
  verboseCmd->SetParameterName("level",false);
}

G4SDmessenger::~G4SDmessenger()
{
  delete listCmd;
  delete activeCmd;
  delete inactiveCmd;
  delete verboseCmd;
  delete hitsDir;
}

void G4SDmessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if( command==listCmd )
  { fSDMan->ListTree(); }
  if( command==activeCmd )
  { fSDMan->Activate(newVal,1); }
  if( command==inactiveCmd )
  { fSDMan->Activate(newVal,0); }
  if( command==verboseCmd )
  { fSDMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal)); }
  return;
}


