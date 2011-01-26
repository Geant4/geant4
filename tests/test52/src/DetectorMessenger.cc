
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


DetectorMessenger::DetectorMessenger(DetectorConstruction* detConstr) :
    detectorConstruction(detConstr) {

  targetDirectory = new G4UIdirectory("/target/");
  targetDirectory -> SetGuidance("Target related commands");

  frontLayerCmd = 
           new G4UIcmdWithAString("/target/frontLayer",this);
  frontLayerCmd -> SetGuidance("Specification of additional front layer");
  frontLayerCmd -> SetParameterName("frontLayer", false);
  frontLayerCmd -> AvailableForStates(G4State_Idle);

  layerRadiusCmd = 
           new G4UIcmdWithADoubleAndUnit("/target/radius",this);
  layerRadiusCmd -> SetGuidance("Specification of layer radius");
  layerRadiusCmd -> SetParameterName("layerRadius", false);
  layerRadiusCmd -> SetRange("layerRadius>0.");
  layerRadiusCmd -> SetUnitCategory("Length");
  layerRadiusCmd -> AvailableForStates(G4State_Idle);

  layerThicknessCmd = 
           new G4UIcmdWithADoubleAndUnit("/target/layerThickness",this);
  layerThicknessCmd -> SetGuidance("Specification of layer thickness");
  layerThicknessCmd -> SetParameterName("layerThickness", false);
  layerThicknessCmd -> SetRange("layerThickness>0.");
  layerThicknessCmd -> SetUnitCategory("Length");
  layerThicknessCmd -> AvailableForStates(G4State_Idle);

  layerMaterialCmd = 
           new G4UIcmdWithAString("/target/layerMaterial",this);
  layerMaterialCmd -> SetGuidance("Specification of layer material");
  layerMaterialCmd -> SetParameterName("layerMaterial", false);
  layerMaterialCmd -> AvailableForStates(G4State_Idle);

  layerMaxStepSizeCmd = 
           new G4UIcmdWithADoubleAndUnit("/target/maxStepSize",this);
  layerMaxStepSizeCmd -> SetGuidance("Specific. of max. step size in target");
  layerMaxStepSizeCmd -> SetParameterName("maxStepSize", false);
  layerMaxStepSizeCmd -> SetRange("maxStepSize>0.");
  layerMaxStepSizeCmd -> SetUnitCategory("Length");
  layerMaxStepSizeCmd -> AvailableForStates(G4State_Idle);

  calPositionCmd = 
           new G4UIcmdWithADoubleAndUnit("/target/calPosition",this);
  calPositionCmd -> SetGuidance("Specific. of a calorimeter at given position");
  calPositionCmd -> SetParameterName("calPosition", false);
  calPositionCmd -> SetUnitCategory("Length");
  calPositionCmd -> AvailableForStates(G4State_Idle);

  calThicknessCmd = 
           new G4UIcmdWithADoubleAndUnit("/target/calThickness",this);
  calThicknessCmd -> SetGuidance("Specific. of a calorimeter at given position");
  calThicknessCmd -> SetParameterName("calThickness", false);
  calThicknessCmd -> SetUnitCategory("Length");
  calThicknessCmd -> AvailableForStates(G4State_Idle);
}


DetectorMessenger::~DetectorMessenger() {
 
  delete calThicknessCmd;
  delete calPositionCmd;
  delete layerMaxStepSizeCmd;
  delete layerMaterialCmd;
  delete layerThicknessCmd;
  delete layerRadiusCmd;
  delete frontLayerCmd;
  delete targetDirectory;
}


void DetectorMessenger::SetNewValue(G4UIcommand* cmd, G4String val) {

  if(cmd == frontLayerCmd) {
     detectorConstruction -> CreateFrontLayer(val);
  }
  if(cmd == layerRadiusCmd) {
     detectorConstruction -> SetLayerRadius(
                               layerRadiusCmd -> GetNewDoubleValue(val));
  }
  if(cmd == layerThicknessCmd) {
     detectorConstruction -> SetLayerThickness(
                               layerThicknessCmd -> GetNewDoubleValue(val));
  } 
  if(cmd == layerMaterialCmd) {
     detectorConstruction -> SetLayerMaterial(val);
  } 
  if(cmd == layerMaxStepSizeCmd) {
     detectorConstruction -> SetLayerMaxStepSize(
                               layerMaxStepSizeCmd -> GetNewDoubleValue(val));
  }
  if(cmd == calPositionCmd) {
     detectorConstruction -> CreateCalorimeter(
                               calPositionCmd -> GetNewDoubleValue(val));
  }
  if(cmd == calThicknessCmd) {
     detectorConstruction -> SetCalorimeterThickness(
                               calThicknessCmd -> GetNewDoubleValue(val));
  }

}
