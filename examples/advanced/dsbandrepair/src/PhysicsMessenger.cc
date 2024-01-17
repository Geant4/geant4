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
/// \file PhysicsMessenger.cc
/// \brief Implementation of the PhysicsMessenger class

#include "PhysicsMessenger.hh"
#include "PhysicsList.hh"

PhysicsMessenger::PhysicsMessenger(PhysicsList* phys)
: G4UImessenger(), fPhysicsList(phys)
{
    fPhysDir = std::make_unique<G4UIdirectory>("/dsbandrepair/phys/");
    fPhysDir->SetGuidance("physics list commands");

    fPhysListCmd = std::make_unique<G4UIcmdWithAString>("/dsbandrepair/phys/physicsList",this);
    fPhysListCmd->SetGuidance("Add modula physics list.");
    fPhysListCmd->SetParameterName("physList",false);
    fPhysListCmd->AvailableForStates(G4State_PreInit);
    fPhysListCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsMessenger::SetNewValue(G4UIcommand* cmd,G4String newValue)
{
    if (cmd == fPhysListCmd.get()) {
        fPhysicsList->RegisterPhysicsList(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......