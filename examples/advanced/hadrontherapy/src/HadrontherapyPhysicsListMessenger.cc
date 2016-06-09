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
// $Id: HadrontherapyPhisicsListMessenger.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
#include "HadrontherapyPhysicsListMessenger.hh"
#include "HadrontherapyPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

HadrontherapyPhysicsListMessenger::HadrontherapyPhysicsListMessenger(HadrontherapyPhysicsList * physList)
:physicsList(physList)
{  
 listDir = new G4UIdirectory("/physics/");
  // Building modular PhysicsList

 physicsListCmd = new G4UIcmdWithAString("/physics/addPhysics",this);  
 physicsListCmd->SetGuidance("Add chunks of PhysicsList.");
 physicsListCmd->SetParameterName("physList",false);
 physicsListCmd->AvailableForStates(G4State_PreInit);
}

HadrontherapyPhysicsListMessenger::~HadrontherapyPhysicsListMessenger()
{
  delete physicsListCmd;
  delete listDir;
}

void HadrontherapyPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if (command == physicsListCmd)
    { physicsList->AddPhysicsList(newValue);} 
}






