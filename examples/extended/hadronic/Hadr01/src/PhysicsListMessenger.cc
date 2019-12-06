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
/// \file hadronic/Hadr01/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// PhysicsListMessenger
//
// Created: 31.01.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of Hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//
// 

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:G4UImessenger(), fPhysicsList(pPhys)
{   
  fPListCmd = new G4UIcmdWithAString("/testhadr/Physics",this);
  fPListCmd->SetGuidance("Add modular physics list.");
  fPListCmd->SetParameterName("PList",false);
  fPListCmd->AvailableForStates(G4State_PreInit);

  fListCmd = new G4UIcmdWithoutParameter("/testhadr/ListPhysics",this);
  fListCmd->SetGuidance("Available Physics Lists");
  fListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fPListCmd;
  delete fListCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fPListCmd ) {
    if(fPhysicsList) {
      G4String name = newValue;
      if(name == "PHYSLIST") {
        char* path = std::getenv(name);
        if (path) name = G4String(path);
        else {
          G4cout << "### PhysicsListMessenger WARNING: "
                 << " environment variable PHYSLIST is not defined"
                 << G4endl;
          return; 
        }
      }
      fPhysicsList->AddPhysicsList(name);
    } else {
      G4cout << "### PhysicsListMessenger WARNING: "
             << " /testhadr/Physics UI command is not available "
             << "for reference Physics List" << G4endl;
    }

  } else if( command == fListCmd ) {
    if(fPhysicsList) {
      fPhysicsList->List();
    } else { 
      G4cout << "### PhysicsListMessenger WARNING: "
             << " /testhadr/ListPhysics UI command is not available "
             << "for reference Physics List" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
