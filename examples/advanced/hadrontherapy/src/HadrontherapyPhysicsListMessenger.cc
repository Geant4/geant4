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
// $Id: HadrontherapyPhisicsListMessenger.cc; Nov 2008
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a), M.P. Russo
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

 packageListCmd = new G4UIcmdWithAString("/physics/addPackage",this);
 packageListCmd->SetGuidance("Add physics package.");
 packageListCmd->SetParameterName("package",false);
 packageListCmd->AvailableForStates(G4State_PreInit);
}

HadrontherapyPhysicsListMessenger::~HadrontherapyPhysicsListMessenger()
{
  delete physicsListCmd;
  delete listDir;
  delete packageListCmd;
}

void HadrontherapyPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if (command == physicsListCmd)
   { physicsList->AddPhysicsList(newValue);}
 else if (command == packageListCmd)
   { physicsList->AddPackage(newValue);}
}






