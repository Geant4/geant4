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
// $Id: RemSimPhysicsListMessenger.cc,v 1.3 2004/05/22 12:57:07 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it

#include "RemSimPhysicsListMessenger.hh"
#include "RemSimPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

RemSimPhysicsListMessenger::RemSimPhysicsListMessenger(RemSimPhysicsList * List)
:RemSimList(List)
{
  EnDir = new G4UIdirectory("/physics/");
  EnDir -> SetGuidance("physics commands");

  physicsListCmd = new G4UIcmdWithAString("/physics/addPhysics",this);  
  physicsListCmd -> SetGuidance("Add chunks of PhysicsList");
  physicsListCmd -> SetParameterName("choice",false);
  physicsListCmd -> AvailableForStates(G4State_PreInit);  
}

RemSimPhysicsListMessenger::~RemSimPhysicsListMessenger()
{  
  delete physicsListCmd;
  delete EnDir;
}
  
void RemSimPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == physicsListCmd)
   { RemSimList->AddPhysicsList(newValue); }
}







