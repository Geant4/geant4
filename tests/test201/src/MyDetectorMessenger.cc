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
// $Id: MyDetectorMessenger.cc,v 1.3 2001-07-11 10:10:23 gunter Exp $
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

