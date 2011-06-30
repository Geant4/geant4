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

#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


PhysicsListMessenger::PhysicsListMessenger(PhysicsList* physList) :
    physicsList(physList) {

  physicsDirectory = new G4UIdirectory("/physics/");
  physicsDirectory -> SetGuidance("Physics commands");

  physicsConstrCmd = 
               new G4UIcmdWithAString("/physics/physConstructor", this);  
  physicsConstrCmd -> SetGuidance("Registration of physics constructors");
  physicsConstrCmd -> SetParameterName("physConstructor",false);
  physicsConstrCmd -> AvailableForStates(G4State_PreInit);  

  prodThresholdCmd = 
               new G4UIcmdWithADoubleAndUnit("/physics/prodThreshold", this);
  prodThresholdCmd -> SetGuidance("Specification of production threshold");
  prodThresholdCmd -> SetParameterName("prodThreshold",true);
  prodThresholdCmd -> SetDefaultValue(0.001);
  prodThresholdCmd -> SetDefaultUnit("mm");
  prodThresholdCmd -> AvailableForStates(G4State_PreInit);  
}


PhysicsListMessenger::~PhysicsListMessenger() {

  delete prodThresholdCmd;
  delete physicsConstrCmd;
  delete physicsDirectory;
}


void PhysicsListMessenger::SetNewValue(G4UIcommand* cmd, G4String val) {

  if(cmd == physicsConstrCmd) 
     physicsList -> RegisterPhysConstructor(val);
  if(cmd == prodThresholdCmd) 
     physicsList -> SetProdThreshold(prodThresholdCmd -> GetNewDoubleValue(val));
}
