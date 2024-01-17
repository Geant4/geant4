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
/// \file ChemPhysicsMessenger.cc
/// \brief Implementation of the ChemPhysicsMessenger class

#include "ChemPhysicsMessenger.hh"
#include "ChemPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ChemPhysicsMessenger::ChemPhysicsMessenger(ChemPhysicsList* phys)
: G4UImessenger(), fPhysicsList(phys), fPhysListCmd(nullptr)
{
    fPhysDir = std::make_unique<G4UIdirectory>("/dsbandrepair/chem/");
    fPhysDir->SetGuidance("physics list commands");

    fPhysListCmd = std::make_unique<G4UIcmdWithAString>("/dsbandrepair/chem/physicsList",this);
    fPhysListCmd->SetGuidance("Add modula Physlist.");
    fPhysListCmd->SetParameterName("physList",false);
    fPhysListCmd->AvailableForStates(G4State_PreInit);

    fChemListCmd = std::make_unique<G4UIcmdWithAString>("/dsbandrepair/chem/chemList",this);
    fChemListCmd->SetGuidance("Add modula Chemlist.");
    fChemListCmd->SetParameterName("chemList",false);
    fChemListCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ChemPhysicsMessenger::SetNewValue(G4UIcommand* cmd,G4String newValue)
{
    if (cmd == fPhysListCmd.get()) {
        fPhysicsList->RegisterPhysListConstructor(newValue);
    }

    if (cmd == fChemListCmd.get()) {
        //fPhysicsList->RegisterChemListConstructor(newValue); 
        //L.T.Anh: It seems Chemlist must be put in Physicslist's constructor, 
        //and can not be changed via UI messenger
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......