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
// $Id: WLSPhysicsListMessenger.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/WLSPhysicsListMessenger.cc
/// \brief Implementation of the WLSPhysicsListMessenger class
//
//
#include "globals.hh"

#include "WLSPhysicsListMessenger.hh"
#include "WLSPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PionRadiativeDecayChannel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhysicsListMessenger::WLSPhysicsListMessenger(WLSPhysicsList* pPhys)
  : fPhysicsList(pPhys)
{

    fDirectory = new G4UIdirectory("/WLS/phys/");
    fDirectory->SetGuidance("WLSPhysicsList control");
 
    fSetAbsorptionCMD = new G4UIcmdWithABool("/WLS/setAbsorption", this);
    fSetAbsorptionCMD->SetGuidance("Turn on or off absorption process");
    fSetAbsorptionCMD->AvailableForStates(G4State_Idle);

    fVerboseCmd = new G4UIcmdWithAnInteger("/WLS/phys/verbose",this);
    fVerboseCmd->SetGuidance("set verbose for physics processes");
    fVerboseCmd->SetParameterName("verbose",true);
    fVerboseCmd->SetDefaultValue(1);
    fVerboseCmd->SetRange("verbose>=0");
    fVerboseCmd->AvailableForStates(G4State_Idle);
 
    fCerenkovCmd =
                new G4UIcmdWithAnInteger("/WLS/phys/cerenkovMaxPhotons",this);
    fCerenkovCmd->SetGuidance("set max nb of photons per step");
    fCerenkovCmd->SetParameterName("MaxNumber",false);
    fCerenkovCmd->SetRange("MaxNumber>=0");
    fCerenkovCmd->AvailableForStates(G4State_Idle);

    fGammaCutCMD = new G4UIcmdWithADoubleAndUnit("/WLS/phys/gammaCut",this);
    fGammaCutCMD->SetGuidance("Set gamma cut");
    fGammaCutCMD->SetParameterName("Gcut",false);
    fGammaCutCMD->SetUnitCategory("Length");
    fGammaCutCMD->SetRange("Gcut>0.0");
    fGammaCutCMD->SetDefaultUnit("mm");
    fGammaCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fElectCutCMD = new G4UIcmdWithADoubleAndUnit("/WLS/phys/electronCut",this);
    fElectCutCMD->SetGuidance("Set electron cut");
    fElectCutCMD->SetParameterName("Ecut",false);
    fElectCutCMD->SetUnitCategory("Length");
    fElectCutCMD->SetRange("Ecut>0.0");
    fElectCutCMD->SetDefaultUnit("mm");
    fElectCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fPosCutCMD = new G4UIcmdWithADoubleAndUnit("/WLS/phys/positronCut",this);
    fPosCutCMD->SetGuidance("Set positron cut");
    fPosCutCMD->SetParameterName("Pcut",false);
    fPosCutCMD->SetUnitCategory("Length");
    fPosCutCMD->SetRange("Pcut>0.0");
    fPosCutCMD->SetDefaultUnit("mm");
    fPosCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fAllCutCMD = new G4UIcmdWithADoubleAndUnit("/WLS/phys/allCuts",this);
    fAllCutCMD->SetGuidance("Set cut for all");
    fAllCutCMD->SetParameterName("cut",false);
    fAllCutCMD->SetUnitCategory("Length");
    fAllCutCMD->SetRange("cut>0.0");
    fAllCutCMD->SetDefaultUnit("mm");
    fAllCutCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fStepMaxCMD = new G4UIcmdWithADoubleAndUnit("/WLS/phys/stepMax",this);
    fStepMaxCMD->SetGuidance("Set max. step length in the detector");
    fStepMaxCMD->SetParameterName("mxStep",false);
    fStepMaxCMD->SetUnitCategory("Length");
    fStepMaxCMD->SetRange("mxStep>0.0");
    fStepMaxCMD->SetDefaultUnit("mm");
    fStepMaxCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fClearPhysicsCMD =
                  new G4UIcmdWithoutParameter("/WLS/phys/clearPhysics",this);
    fClearPhysicsCMD->SetGuidance("Clear the physics list");
    fClearPhysicsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fRemovePhysicsCMD = new G4UIcmdWithAString("/WLS/phys/removePhysics",this);
    fRemovePhysicsCMD->
                     SetGuidance("Remove a physics process from Physics List");
    fRemovePhysicsCMD->SetParameterName("PList",false);
    fRemovePhysicsCMD->AvailableForStates(G4State_PreInit,G4State_Idle);

    fListCMD = new G4UIcmdWithoutParameter("/WLS/phys/list",this);
    fListCMD->SetGuidance("Available Physics Lists");
    fListCMD->AvailableForStates(G4State_Idle);

    fDecayDirectory = new G4UIdirectory("/decay/");
    fDecayDirectory->SetGuidance("Decay chain control commands.");

    fPienuCMD = new G4UIcmdWithoutParameter("/decay/pienu", this);
    fPienuCMD->SetGuidance("Sets the pi+ to decay into e+, nu");

    fPimunuCMD = new G4UIcmdWithoutParameter("/decay/pimunu", this);
    fPimunuCMD->SetGuidance("Sets the pi+ to decay into mu+, nu");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhysicsListMessenger::~WLSPhysicsListMessenger()
{
    delete fVerboseCmd;
    delete fCerenkovCmd;

    delete fSetAbsorptionCMD;

    delete fGammaCutCMD;
    delete fElectCutCMD;
    delete fPosCutCMD;
    delete fAllCutCMD;

    delete fClearPhysicsCMD;
    delete fRemovePhysicsCMD;

    delete fListCMD;

    delete fPienuCMD;
    delete fPimunuCMD;

    delete fDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
    if( command == fSetAbsorptionCMD ) {
       fPhysicsList->SetAbsorption(G4UIcmdWithABool::GetNewBoolValue(newValue));
    }

    else if( command == fVerboseCmd ) {
       fPhysicsList->SetVerbose(fVerboseCmd->GetNewIntValue(newValue));
    }

    else if( command == fCerenkovCmd ) {
       fPhysicsList->
           SetNbOfPhotonsCerenkov(fCerenkovCmd->GetNewIntValue(newValue));
    }

    else if (command == fPienuCMD) {
       G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
       G4ParticleDefinition* particleDef = particleTable->FindParticle("pi+");
       G4VDecayChannel* mode = 
                        new G4PhaseSpaceDecayChannel("pi+",1.0,2,"e+","nu_e");
       G4DecayTable* table = new G4DecayTable();
       table->Insert(mode);
      // mode = new G4PionRadiativeDecayChannel("pi+",0.000017);
      // table->Insert(mode);
       particleDef->SetDecayTable(table);
    }

    else if (command == fPimunuCMD) {
       G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
       G4ParticleDefinition* particleDef = particleTable->FindParticle("pi+");
       G4VDecayChannel* mode =
                     new G4PhaseSpaceDecayChannel("pi+",1.000,2,"mu+","nu_mu");
       G4DecayTable* table = new G4DecayTable();
       table->Insert(mode);
       particleDef->SetDecayTable(table);
    }

    else if (command == fGammaCutCMD) {
        fPhysicsList->SetCutForGamma(fGammaCutCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fElectCutCMD) {
        fPhysicsList->SetCutForElectron(fElectCutCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fPosCutCMD) {
        fPhysicsList->SetCutForPositron(fPosCutCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fAllCutCMD) {
        G4double cut = fAllCutCMD->GetNewDoubleValue(newValue);
        fPhysicsList->SetCutForGamma(cut);
        fPhysicsList->SetCutForElectron(cut);
        fPhysicsList->SetCutForPositron(cut);
    }
    else if (command == fStepMaxCMD) {
        fPhysicsList->SetStepMax(fStepMaxCMD
                                     ->GetNewDoubleValue(newValue));
    }
    else if (command == fClearPhysicsCMD) {
        fPhysicsList->ClearPhysics();
    }
    else if (command == fRemovePhysicsCMD) {
        G4String name = newValue;
        fPhysicsList->RemoveFromPhysicsList(name);
    }
}
