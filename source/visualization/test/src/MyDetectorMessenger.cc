// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorMessenger.cc,v 1.1 1999-04-16 10:32:36 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyDetectorMessenger.hh"

#include "MyDetectorConstruction.hh"
#include "G4UIparameter.hh"
#include "globals.hh"

MyDetectorMessenger::MyDetectorMessenger(MyDetectorConstruction * myDet)
:myDetector(myDet)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand("/mydet/",this);
  command->SetGuidance("My detector control.");
  fpMyDetCommandDirectory = command;

  command = new G4UIcommand("/mydet/calMaterial",this);
  command->SetGuidance("Define material of the calorimeter.");
  command->SetGuidance(" Available materials :");
  command->SetGuidance("   Air (defaust), Al, Fe, Pb");
  param = new G4UIparameter("material",'c',true);
  param->SetDefaultValue("Air");
  command->SetParameter(param);
  fpCalMaterialCommand = command;
}

MyDetectorMessenger::~MyDetectorMessenger () {
  delete fpCalMaterialCommand;
  delete fpMyDetCommandDirectory;
}

void MyDetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if (command == fpCalMaterialCommand)
  {
    myDetector->SetCalMaterial(newValues);
  }
}

