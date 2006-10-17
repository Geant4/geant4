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
//    *****************************************
//    *                                       *
//    *        CellPhysicsListMessenger.cc    *
//    *                                       *
//    *****************************************
//
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "CellPhysicsListMessenger.hh"
#include "CellPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

CellPhysicsListMessenger::CellPhysicsListMessenger(CellPhysicsList * List)
:CellList(List)
{
  EnDir = new G4UIdirectory("/physics/");
  EnDir -> SetGuidance("physics commands");

  // Create the interactive command to change the cut
  cutECmd = new G4UIcmdWithADoubleAndUnit("/physics/cutE",this);
  cutECmd -> SetGuidance("Set cut values.");
  cutECmd -> SetParameterName("range",true);
  cutECmd -> SetDefaultValue(0.1);
  cutECmd -> SetDefaultUnit("mm");
  cutECmd -> AvailableForStates(G4State_PreInit);
}

CellPhysicsListMessenger::~CellPhysicsListMessenger()
{  
  delete cutECmd;
  delete EnDir;
}
  
void CellPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  // If the cut is changed interactively in idle session
  // this change is notified to the physics component
 if (command == cutECmd)
   {CellList -> SetParticleCut(cutECmd->GetNewDoubleValue(newValue)); }
}







