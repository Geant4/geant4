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

#include "ExExChDetectorConstructionMessenger.hh"
#include "ExExChDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

ExExChDetectorConstructionMessenger::
ExExChDetectorConstructionMessenger(
        ExExChDetectorConstruction* mpga)
:fTarget(mpga){
    fMyXtalDirectory = new G4UIdirectory("/xtal/");
    fMyXtalDirectory->SetGuidance("Xtal setup control commands.");

    fXtalMaterialCmd = new G4UIcmdWithAString("/xtal/setMaterial",
                                              this);
    fXtalMaterialCmd->SetGuidance("Set Xtal material.");
    fXtalMaterialCmd->SetParameterName("xMat",true);
    fXtalMaterialCmd->SetDefaultValue("G4_Si");
    
    fXtalCurvatureRadiusCmd =
        new G4UIcmdWith3VectorAndUnit("/xtal/setCurvRadius",this);
    fXtalCurvatureRadiusCmd->SetGuidance("Set Xtal curvature radius.");
    fXtalCurvatureRadiusCmd->SetParameterName("xtalcrx",
                                              "xtalcry",
                                              "xtalcrz",
                                              true);
    fXtalCurvatureRadiusCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalCurvatureRadiusCmd->SetDefaultUnit("m");

    fXtalSizeCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setSize",this);
    fXtalSizeCmd->SetGuidance("Set Xtal size.");
    fXtalSizeCmd->SetParameterName("xtalSizeX",
                                   "xtalSizeY",
                                   "xtalSizeZ",
                                   true);
    fXtalSizeCmd->SetDefaultValue(G4ThreeVector(6.,2.,6.));
    fXtalSizeCmd->SetDefaultUnit("mm");
    
    fXtalAngleCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setAngle",this);
    fXtalAngleCmd->SetGuidance("Set Xtal angles in the world.");
    fXtalAngleCmd->SetParameterName("xtalAngleX",
                                    "xtalAngleY",
                                    "xtalAngleZ",
                                    true);
    fXtalAngleCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
    fXtalAngleCmd->SetDefaultUnit("rad");
    
    fXtalCellSizeCmd = new G4UIcmdWith3VectorAndUnit("/xtal/setCellSize",this);
    fXtalCellSizeCmd->SetGuidance("Set Xtal cell size.");
    fXtalCellSizeCmd->SetParameterName("xtalCellSizeX",
                                       "xtalCellSizeY",
                                       "xtalCellSizeZ",
                                       true);
    fXtalCellSizeCmd->SetDefaultValue(G4ThreeVector(5.431,5.431,5.431));
    fXtalCellSizeCmd->SetDefaultUnit("angstrom");
    
    fXtalCellAngleCmd =
        new G4UIcmdWith3VectorAndUnit("/xtal/setCellAngle",this);
    fXtalCellAngleCmd->SetGuidance("Set Xtal cell angles in the world.");
    fXtalCellAngleCmd->SetParameterName("xtalCellAngleX",
                                        "xtalCellAngleY",
                                        "xtalCellAngleZ",
                                        true);
    fXtalCellAngleCmd->SetDefaultValue(G4ThreeVector(90.,90.,90.));
    fXtalCellAngleCmd->SetDefaultUnit("deg");

    fXtalCellThermalVibration =
        new G4UIcmdWithADoubleAndUnit("/xtal/setThermVibr",this);
    fXtalCellThermalVibration->SetGuidance("Thermal vibration amplitude");
    fXtalCellThermalVibration->SetParameterName("xtalThermVibr",true);
    fXtalCellThermalVibration->SetDefaultValue(0.075);
    fXtalCellThermalVibration->SetDefaultUnit("angstrom");

    fXtalMillerCmd = new G4UIcmdWith3Vector("/xtal/setMiller",this);
    fXtalMillerCmd->SetGuidance("Set Miller Indexes");
    fXtalMillerCmd->SetParameterName("xtalMiller1",
                                     "xtalMiller2",
                                     "xtalMiller3",
                                     true);
    fXtalMillerCmd->SetDefaultValue(G4ThreeVector(2,2,0));
    
    fAddXtal = new G4UIcmdWithAString("/addxtal",this);
    fAddXtal->SetGuidance("Add Xtal.");
    fAddXtal->SetParameterName("addxtal",true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExExChDetectorConstructionMessenger::
~ExExChDetectorConstructionMessenger(){
    delete fXtalMaterialCmd;
    delete fXtalMillerCmd;
    delete fXtalCurvatureRadiusCmd;
    delete fXtalSizeCmd;
    delete fXtalAngleCmd;
    delete fXtalCellSizeCmd;
    delete fXtalCellAngleCmd;
    delete fAddXtal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExExChDetectorConstructionMessenger::SetNewValue(
                                                    G4UIcommand *command,
                                                    G4String newValue){
    if(command==fXtalMaterialCmd ){
        fTarget->SetXtalMaterial(newValue);
    }
    if(command==fXtalCurvatureRadiusCmd ){
        fTarget->SetXtalCurvatureRadius(
                        fXtalCurvatureRadiusCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalSizeCmd ){
        fTarget->SetXtalSize(fXtalSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalAngleCmd ){
        fTarget->SetXtalAngle(fXtalAngleCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellSizeCmd ){
        fTarget->SetXtalCellSize(
                        fXtalCellSizeCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellAngleCmd ){
        fTarget->SetXtalCellAngle(
                        fXtalCellAngleCmd->GetNew3VectorValue(newValue));
    }
    if(command==fXtalCellThermalVibration ){
        fTarget->SetXtalThermalVibrationAmplitude(
                        fXtalCellThermalVibration->GetNewDoubleValue(newValue));
    }
    if(command==fXtalMillerCmd ){
        fTarget->SetXtalMiller(fXtalMillerCmd->GetNew3VectorValue(newValue));
    }
    if(command==fAddXtal ){
        fTarget->AddXtalTarget();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String ExExChDetectorConstructionMessenger::GetCurrentValue(
                                                    G4UIcommand * command){
    G4String cv;
    
    if( command==fXtalMaterialCmd ){
        cv = fTarget->GetXtalMaterial();
    }
    if( command==fXtalCurvatureRadiusCmd ){
        cv = fXtalCurvatureRadiusCmd->ConvertToString(
                fTarget->GetXtalCurvatureRadius(),"m");
    }
    if( command==fXtalSizeCmd ){
        cv = fXtalSizeCmd->ConvertToString(fTarget->GetXtalSize(),"mm");
    }
    if( command==fXtalAngleCmd ){
        cv = fXtalAngleCmd->ConvertToString(fTarget->GetXtalAngle(),"rad");
    }
    if( command==fXtalCellSizeCmd ){
        cv = fXtalCellSizeCmd->ConvertToString(
                fTarget->GetXtalCellSize(),"angstrom");
    }
    if( command==fXtalCellAngleCmd ){
        cv = fXtalCellAngleCmd->ConvertToString(
                fTarget->GetXtalCellAngle(),"deg");
    }
    if( command==fXtalCellThermalVibration ){
        cv = fXtalCellThermalVibration->ConvertToString(
                fTarget->GetXtalThermalVibrationAmplitude(),"angstrom");
    }
    if( command==fXtalMillerCmd ){
        cv = fXtalMillerCmd->ConvertToString(fTarget->GetXtalMiller());
    }
    
    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
