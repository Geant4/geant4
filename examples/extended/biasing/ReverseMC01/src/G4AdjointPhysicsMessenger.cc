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
/// \file biasing/ReverseMC01/src/G4AdjointPhysicsMessenger.cc
/// \brief Implementation of the G4AdjointPhysicsMessenger class
//
//
//////////////////////////////////////////////////////////////
//      Class Name:        G4AdjointPhysicsMessenger
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AdjointPhysicsMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UnitsTable.hh"
#include "G4AdjointPhysicsList.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointPhysicsMessenger::G4AdjointPhysicsMessenger(
                                           G4AdjointPhysicsList* pPhysicsList)
: G4UImessenger(),
  fPhysicsList(pPhysicsList),
  fPhysicsDir(0),
  fUsepIonisationCmd(0),
  fUseBremCmd(0),
  fUseComptonCmd(0),
  fUseMSCmd(0),
  fUsePEEffectCmd(0),
  fUseGammaConversionCmd(0),
  fUseEgainFluctuationCmd(0),
  fSetEminAdjModelsCmd(0),
  fSetEmaxAdjModelsCmd(0)
{ 
  fPhysicsDir = new G4UIdirectory("/adjoint_physics/");

  fPhysicsDir->SetGuidance(
                     "Definition of adjoint and forward physics processes");
  //-------
  fUsepIonisationCmd = new G4UIcmdWithABool(
                 "/adjoint_physics/UseProtonIonisation",this);
  fUsepIonisationCmd->SetGuidance(
      "If true (false) the proton ionisation is (not) considered");
  fUsepIonisationCmd->AvailableForStates(G4State_PreInit);
  
  fUseBremCmd = new G4UIcmdWithABool("/adjoint_physics/UseBremsstrahlung",this);
  fUseBremCmd->SetGuidance(
          "If true (false) the bremsstrahlung process is (not) considered");
  fUseBremCmd->AvailableForStates(G4State_PreInit); 
  
  fUseComptonCmd = new G4UIcmdWithABool("/adjoint_physics/UseCompton",this);
  fUseComptonCmd->SetGuidance(
                "If true (false) the Compton scattering is (not) considered");
  fUseComptonCmd->AvailableForStates(G4State_PreInit);  
  
  fUseMSCmd = new G4UIcmdWithABool("/adjoint_physics/UseMS",this);
  fUseMSCmd->SetGuidance(
      "If true (false) the continuous multiple scattering is (not) considered");
  fUseMSCmd->AvailableForStates(G4State_PreInit); 
  
  fUseEgainFluctuationCmd = new G4UIcmdWithABool(
              "/adjoint_physics/UseEgainElossFluctuation",this);
  fUseEgainFluctuationCmd->SetGuidance(
    "Switch on/off the fluctation for continuous energy gain/loss");
  fUseEgainFluctuationCmd->AvailableForStates(G4State_PreInit); 
  
  fUsePEEffectCmd = new G4UIcmdWithABool("/adjoint_physics/UsePEEffect",this);
  fUsePEEffectCmd->AvailableForStates(G4State_PreInit); 
  fUsePEEffectCmd->SetGuidance(
     "If true (false) the photo electric effect is (not) considered");
   
  fUseGammaConversionCmd = new G4UIcmdWithABool(
          "/adjoint_physics/UseGammaConversion",this);
  fUseGammaConversionCmd->AvailableForStates(G4State_PreInit);
  fUseGammaConversionCmd->SetGuidance(
             "If true the fwd gamma pair conversion is considered");
  
  fSetEminAdjModelsCmd =  new G4UIcmdWithADoubleAndUnit(
             "/adjoint_physics/SetEminForAdjointModels",this);
  fSetEminAdjModelsCmd->SetGuidance(
            "Set the minimum energy  of the adjoint models");
  fSetEminAdjModelsCmd->SetParameterName("Emin",false);
  fSetEminAdjModelsCmd->SetUnitCategory("Energy");
  fSetEminAdjModelsCmd->AvailableForStates(G4State_PreInit);
  
  fSetEmaxAdjModelsCmd =  new G4UIcmdWithADoubleAndUnit(
                 "/adjoint_physics/SetEmaxForAdjointModels",this);
  fSetEmaxAdjModelsCmd->SetGuidance(
                "Set the minimum energy of the adjoint models.");
  fSetEmaxAdjModelsCmd->SetParameterName("Emax",false);
  fSetEmaxAdjModelsCmd->SetUnitCategory("Energy");
  fSetEmaxAdjModelsCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointPhysicsMessenger::~G4AdjointPhysicsMessenger()
{
  delete fUsepIonisationCmd;
  delete fUseBremCmd;
  delete fUseComptonCmd;
  delete fUseMSCmd;
  delete fUsePEEffectCmd;
  delete fUseGammaConversionCmd;
  delete fUseEgainFluctuationCmd;
  delete fSetEminAdjModelsCmd;
  delete fSetEmaxAdjModelsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsMessenger::SetNewValue(G4UIcommand* command,
                                                             G4String newValue)
{
  if ( command==fUsepIonisationCmd){
    fPhysicsList->SetUseProtonIonisation(
                               fUsepIonisationCmd->GetNewBoolValue(newValue));
     
  }
  else if ( command==fUseBremCmd){
          fPhysicsList->SetUseBrem(fUseBremCmd->GetNewBoolValue(newValue));
  }
  else if ( command==fUseComptonCmd){
         fPhysicsList->SetUseCompton(fUseComptonCmd->GetNewBoolValue(newValue));
  }
  else if ( command==fUseMSCmd){
    fPhysicsList->SetUseMS(fUseMSCmd->GetNewBoolValue(newValue));
  }
  else if ( command==fUsePEEffectCmd){
    fPhysicsList->SetUsePEEffect(fUsePEEffectCmd->GetNewBoolValue(newValue));
  }
  else if ( command==fUseGammaConversionCmd){
    fPhysicsList->SetUseGammaConversion(
                         fUseGammaConversionCmd->GetNewBoolValue(newValue));
  }
  else if ( command==fUseEgainFluctuationCmd){
    fPhysicsList->SetUseEgainFluctuation(
                    fUseEgainFluctuationCmd->GetNewBoolValue(newValue));
  }

  else if ( command== fSetEminAdjModelsCmd){
    fPhysicsList->SetEminAdjModels(
                     fSetEminAdjModelsCmd->GetNewDoubleValue(newValue));
  }
  else if ( command== fSetEmaxAdjModelsCmd){
    fPhysicsList->SetEmaxAdjModels(
                      fSetEmaxAdjModelsCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

