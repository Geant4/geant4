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
// $Id: G4SDmessenger.cc 67992 2013-03-13 10:59:57Z gcosmo $
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


