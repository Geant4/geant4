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

#include "ExExChPhysicsListMessenger.hh"
#include "ExExChPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ios.hh"

ExExChPhysicsListMessenger::ExExChPhysicsListMessenger(
                                    ExExChPhysicsList * mpga)
:fTarget (mpga){
    fMyDirectory = new G4UIdirectory("/xtal/");
    
    fFilePotentialNameCmd = new G4UIcmdWithAString("/xtal/potfilename",this);
    fFilePotentialNameCmd->SetGuidance("Filename for input potential vector.");
    fFilePotentialNameCmd->SetParameterName("potfilename",true);
    fFilePotentialNameCmd->SetDefaultValue("");
    fTransverseVariationMaxCmd =
        new G4UIcmdWithADoubleAndUnit("/xtal/setTransVarMax",this);
    fTransverseVariationMaxCmd->SetGuidance("Integration - transverse variation max");
    fTransverseVariationMaxCmd->SetParameterName("trvarmax",true);
    fTransverseVariationMaxCmd->SetDefaultValue(2.E-2);
    fTransverseVariationMaxCmd->SetDefaultUnit("angstrom");
    fTransverseVariationMaxCmd->SetRange("trvarmax>=0.0");


    fTimeStepMinCmd =
        new G4UIcmdWithADoubleAndUnit("/xtal/setTimeStepMin",this);
    fTimeStepMinCmd->SetGuidance("Integration - time step min");
    fTimeStepMinCmd->SetParameterName("timestmin",true);
    fTimeStepMinCmd->SetDefaultValue(2.E0);
    fTimeStepMinCmd->SetDefaultUnit("angstrom");
    fTimeStepMinCmd->SetRange("timestmin>=0.0");}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExExChPhysicsListMessenger::~ExExChPhysicsListMessenger(){
    delete fFilePotentialNameCmd;
    delete fTransverseVariationMaxCmd;
    delete fTimeStepMinCmd;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChPhysicsListMessenger::SetNewValue(G4UIcommand * command,
                                                   G4String newValue){
    if( command==fFilePotentialNameCmd ){
        fTarget->SetFilePotentialName(newValue);
    }
    if(command==fTransverseVariationMaxCmd ){
        fTarget->SetTransverseVariationMax(fTransverseVariationMaxCmd->GetNewDoubleValue(newValue));
    }
    if(command==fTimeStepMinCmd ){
        fTarget->SetTimeStepMin(fTimeStepMinCmd->GetNewDoubleValue(newValue));
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ExExChPhysicsListMessenger::GetCurrentValue(
                                            G4UIcommand * command){
    G4String cv;
    if( command==fFilePotentialNameCmd ){
        cv = fTarget->GetFilePotentialName();
    }
    if(command==fTransverseVariationMaxCmd ){
        cv = fTransverseVariationMaxCmd->ConvertToString(fTarget->GetTransverseVariationMax());
    }
    if(command==fTimeStepMinCmd ){
        cv = fTimeStepMinCmd->ConvertToString(fTarget->GetTimeStepMin());
    }
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
