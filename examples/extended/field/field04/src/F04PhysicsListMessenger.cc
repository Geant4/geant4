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
/// \file field/field04/src/F04PhysicsListMessenger.cc
/// \brief Implementation of the F04PhysicsListMessenger class
//

#include "globals.hh"

#include "F04PhysicsListMessenger.hh"
#include "F04PhysicsList.hh"

#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PionRadiativeDecayChannel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PhysicsListMessenger::F04PhysicsListMessenger(F04PhysicsList* pPhys)
    : fPhysicsList(pPhys)
{
    fDirectory = new G4UIdirectory("/exp/phys/");
    fDirectory->SetGuidance("Control the physics lists");

    fStepMaxCMD = new G4UIcmdWithADoubleAndUnit("/exp/phys/stepMax",this);
    fStepMaxCMD->SetGuidance("Set max. step length in the detector");
    fStepMaxCMD->SetParameterName("mxStep",false);
    fStepMaxCMD->SetUnitCategory("Length");
    fStepMaxCMD->SetRange("mxStep>0.0");
    fStepMaxCMD->SetDefaultUnit("mm");
    fStepMaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
/*
    fClearPhysicsCMD = new G4UIcmdWithoutParameter("/exp/phys/clearPhysics",
                                                                        this);
    fClearPhysicsCMD->SetGuidance("Clear the physics list");
    fClearPhysicsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fRemovePhysicsCMD = new G4UIcmdWithAString("/exp/phys/removePhysics",this);
    fRemovePhysicsCMD->
                     SetGuidance("Remove a physics process from Physics List");
    fRemovePhysicsCMD->SetParameterName("PList",false);
    fRemovePhysicsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);
*/
    fDecayDirectory = new G4UIdirectory("/decay/");
    fDecayDirectory->SetGuidance("Decay chain control commands.");

    fPienuCMD = new G4UIcmdWithoutParameter("/decay/pienu", this);
    fPienuCMD->SetGuidance("Sets the pi+ to decay into e+, nu");

    fPimunuCMD = new G4UIcmdWithoutParameter("/decay/pimunu", this);
    fPimunuCMD->SetGuidance("Sets the pi+ to decay into mu+, nu");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PhysicsListMessenger::~F04PhysicsListMessenger()
{
/*
    delete fClearPhysicsCMD;
    delete fRemovePhysicsCMD;
*/
    delete fPienuCMD;
    delete fPimunuCMD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
    G4ParticleTable* fParticleTable = G4ParticleTable::GetParticleTable();

    if (command == fPienuCMD) {
       G4ParticleDefinition* fParticleDef = fParticleTable->FindParticle("pi+");
       G4VDecayChannel* fMode =
                     new G4PhaseSpaceDecayChannel("pi+",0.999983,2,"e+","nu_e");
       G4DecayTable* fTable = new G4DecayTable();
       fTable->Insert(fMode);
       fMode = new G4PionRadiativeDecayChannel("pi+",0.000017);
       fTable->Insert(fMode);
       fParticleDef->SetDecayTable(fTable);
    }

    if (command == fPimunuCMD) {
       G4ParticleDefinition* fParticleDef = fParticleTable->FindParticle("pi+");
       G4VDecayChannel* fMode =
                     new G4PhaseSpaceDecayChannel("pi+",1.000,2,"mu+","nu_mu");
       G4DecayTable* fTable = new G4DecayTable();
       fTable->Insert(fMode);
       fParticleDef->SetDecayTable(fTable);
    }

    else if (command == fStepMaxCMD) {
        fPhysicsList->SetStepMax(fStepMaxCMD
                                     ->GetNewDoubleValue(newValue));
    }
/*  else if (command == fClearPhysicsCMD) {
        fPhysicsList->ClearPhysics();
    }
    else if (command == fRemovePhysicsCMD) {
        G4String name = newValue;
        fPhysicsList->RemoveFromPhysicsList(name);
    }
*/
}
