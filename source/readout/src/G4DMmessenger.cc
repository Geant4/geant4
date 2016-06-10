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
// $Id: G4DMmessenger.cc 80987 2014-05-19 10:50:22Z gcosmo $
//
// ---------------------------------------------------------------------

#include "G4DMmessenger.hh"
#include "G4DigiManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

G4DMmessenger::G4DMmessenger(G4DigiManager* DigiManager):fDMan(DigiManager)
{
  digiDir = new G4UIdirectory("/digi/");
  digiDir->SetGuidance("DigitizerModule");

  listCmd = new G4UIcmdWithoutParameter("/digi/List",this);
  listCmd->SetGuidance("List names of digitizer modules.");

  digiCmd = new G4UIcmdWithAString("/digi/Digitize",this);
  digiCmd->SetGuidance("Invoke Digitize method of a digitizer module");
  digiCmd->SetParameterName("moduleName",false);

  verboseCmd = new G4UIcmdWithAnInteger("/digi/Verbose",this);
  verboseCmd->SetGuidance("Set the Verbose level.");
  verboseCmd->SetParameterName("level",false);
}

G4DMmessenger::~G4DMmessenger()
{
  delete listCmd;
  delete digiCmd;
  delete verboseCmd;
  delete digiDir;
}

void G4DMmessenger::SetNewValue(G4UIcommand * command,G4String newVal)
{
  if( command==listCmd )
  { fDMan->List(); }
  if( command==digiCmd )
  { fDMan->Digitize(newVal); }
  if( command==verboseCmd )
  { fDMan->SetVerboseLevel(verboseCmd->GetNewIntValue(newVal)); }
  return;
}


