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

#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

DetectorConstructionMessenger::
DetectorConstructionMessenger(
                              DetectorConstruction* mpga)
:fTarget(mpga){
    fMyXtalDirectory = new G4UIdirectory("/xtal/");
    fMyXtalDirectory->SetGuidance("Crystal setup control commands.");
    
    fXtalMaterialCmd = new G4UIcmdWithAString("/xtal/setMaterial",
                                              this);
    fXtalMaterialCmd->SetGuidance("Set crystal material.");
    fXtalMaterialCmd->SetParameterName("xMat",true);
    fXtalMaterialCmd->SetDefaultValue("G4_Si");
    
    
    fXtalSizeCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setSize",this);
    fXtalSizeCmd->SetGuidance("Set crystal size.");
    fXtalSizeCmd->SetParameterName("xtalSizeX",
                                   "xtalSizeY",
                                   "xtalSizeZ",
                                   true);
    fXtalSizeCmd->SetDefaultValue(G4ThreeVector(6.,2.,6.));
    fXtalSizeCmd->SetDefaultUnit("mm");
    
    fXtalBRCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setBR",this);
    fXtalBRCmd->SetGuidance("Set crystal bending radius.");
    fXtalBRCmd->SetParameterName("xtalBRX",
                                 "xtalBRY",
                                 "xtalBRZ",
                                 true);
    fXtalBRCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalBRCmd->SetDefaultUnit("m");
    
    
    fXtalAngleCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setAngle",this);
    fXtalAngleCmd->SetGuidance("Set crystal angles.");
    fXtalAngleCmd->SetParameterName("xtalAngleX",
                                    "xtalAngleY",
                                    "xtalAngleZ",
                                    true);
    fXtalAngleCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalAngleCmd->SetDefaultUnit("rad");
    
    fXtalECCmd = new G4UIcmdWithAString("/xtal/setEC",
                                        this);
    fXtalECCmd->SetGuidance("Set crystal EC.");
    fXtalECCmd->SetParameterName("xEC",true);
    fXtalECCmd->SetDefaultValue("data/Si220pl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::
~DetectorConstructionMessenger(){
    delete fXtalMaterialCmd;
    delete fXtalSizeCmd;
    delete fXtalAngleCmd;
    delete fXtalECCmd;
    delete fXtalBRCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(
                                                G4UIcommand *command,
                                                G4String newValue){
    if(command==fXtalMaterialCmd ){
        fTarget->SetMaterial(newValue);
    }
    if(command==fXtalSizeCmd ){
        fTarget->SetSizes(fXtalSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalBRCmd ){
        fTarget->SetBR(fXtalBRCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalAngleCmd ){
        fTarget->SetAngles(fXtalAngleCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalECCmd ){
        fTarget->SetEC(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String DetectorConstructionMessenger::GetCurrentValue(
                                                        G4UIcommand * command){
    G4String cv;
    
    if( command==fXtalMaterialCmd ){
        cv = fTarget->GetMaterial();
    }
    if( command==fXtalSizeCmd ){
        cv = fXtalSizeCmd->ConvertToString(fTarget->GetSizes(),"mm");
    }
    if( command==fXtalBRCmd ){
        cv = fXtalBRCmd->ConvertToString(fTarget->GetBR(),"m");
    }
    if( command==fXtalAngleCmd ){
        cv = fXtalAngleCmd->ConvertToString(fTarget->GetAngles(),"rad");
    }
    if( command==fXtalECCmd ){
        cv = fTarget->GetEC();
    }
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
