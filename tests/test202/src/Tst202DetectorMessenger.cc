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
// $Id: Tst202DetectorMessenger.cc,v 1.1 2007-02-08 15:45:29 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst202DetectorMessenger.hh"

#include "Tst202DetectorConstruction.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"

Tst202DetectorMessenger::Tst202DetectorMessenger(Tst202DetectorConstruction * myDet)
:myDetector(myDet)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand("/test202/",this);
  command->SetGuidance("Tst202 detector control.");
  fpTst202DetCommandDirectory = command;

  command = new G4UIcommand("/test202/calMaterial",this);
  command->SetGuidance("Define material of the calorimeter.");
  command->SetGuidance(" Available materials :");
  command->SetGuidance("   Air (defaust), Al, Fe, Pb");
  param = new G4UIparameter("material",'c',true);
  param->SetDefaultValue("Air");
  command->SetParameter(param);
  fpCalMaterialCommand = command;

  command = new G4UIcommand("/test202/select",this);
  command->SetGuidance("Select volume");
  param = new G4UIparameter("volume",'s',true);
  param->SetDefaultValue("none");
  command->SetParameter(param);
  param = new G4UIparameter("select",'b',true);
  param->SetDefaultValue(true);
  command->SetParameter(param);
  fpVolumeSelectionCommand = command;

  fpVerboseCommand = new G4UIcmdWithAnInteger("/test202/verbose",this);
  fpVerboseCommand->SetGuidance
    ("Verbosity: 0: silent, >1: various degrees of verbosity");
  fpVerboseCommand->SetParameterName("verbosity",true,false);
  fpVerboseCommand->SetDefaultValue(0);

}

Tst202DetectorMessenger::~Tst202DetectorMessenger () {
  delete fpVolumeSelectionCommand;
  delete fpCalMaterialCommand;
  delete fpTst202DetCommandDirectory;
}

void Tst202DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if (command == fpCalMaterialCommand)
  {
    myDetector->SetCalMaterial(newValues);
  }
  else if (command == fpVolumeSelectionCommand)
  {
    G4String volumeName, selectString;
    std::istringstream iss(newValues);
    iss >> volumeName >> selectString;
    G4bool select = fpVolumeSelectionCommand->ConvertToBool(selectString);
    myDetector->SetVolume(volumeName,select);
  }
  else if (command == fpVerboseCommand)
  {
    myDetector->SetVerbosity(fpVerboseCommand->GetNewIntValue(newValues));
  }
}

