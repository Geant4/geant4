//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorMessenger::XrayFluoDetectorMessenger(XrayFluoDetectorConstruction * Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");
  
  SamMaterCmd = new G4UIcmdWithAString("/apparate/setSamMat",this);
  SamMaterCmd->SetGuidance("Select Material of the Sample.");
  SamMaterCmd->SetParameterName("choice",false);
  SamMaterCmd->AvailableForStates(Idle);

  SenMaterCmd = new G4UIcmdWithAString("/apparate/setSenMat",this);
  SenMaterCmd->SetGuidance("Select Material of the Sensor.");
  SenMaterCmd->SetParameterName("choice",false);
  SenMaterCmd->AvailableForStates(Idle);

  SamThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSamThick",this);
  SamThickCmd->SetGuidance("Set Thickness of the Sample");
  SamThickCmd->SetParameterName("Size",false);
  SamThickCmd->SetRange("Size>=0.");
  SamThickCmd->SetUnitCategory("Length");
  SamThickCmd->AvailableForStates(Idle);
  
  SenThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSenThick",this);
  SenThickCmd->SetGuidance("Set Thickness of the Sensor");
  SenThickCmd->SetParameterName("Size",false);
  SenThickCmd->SetRange("Size>=0.");
  SenThickCmd->SetUnitCategory("Length");  
  SenThickCmd->AvailableForStates(Idle);
 
  SamSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSamSizeYZ",this);
  SamSizeYZCmd->SetGuidance("Set tranverse size of the sample");
  SamSizeYZCmd->SetParameterName("Size",false);
  SamSizeYZCmd->SetRange("Size>0.");
  SamSizeYZCmd->SetUnitCategory("Length");    
  SamSizeYZCmd->AvailableForStates(Idle);

  SenSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSenSizeYZ",this);
  SenSizeYZCmd->SetGuidance("Set tranverse size of the sensor");
  SenSizeYZCmd->SetParameterName("Size",false);
  SenSizeYZCmd->SetRange("Size>0.");
  SenSizeYZCmd->SetUnitCategory("Length");    
  SenSizeYZCmd->AvailableForStates(Idle);

  SenDistCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSenDist",this);
  SenDistCmd->SetGuidance("Set distance of the detector");
  SenDistCmd->SetParameterName("Size",false);
  SenDistCmd->SetRange("Size>0.");
  SenDistCmd->SetUnitCategory("Length");    
  SenDistCmd->AvailableForStates(Idle);

  SenAngleCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSenAngle",this);
  SenAngleCmd->SetGuidance("Set angular displacement of the sensor");
  SenAngleCmd->SetParameterName("Size",false);
  SenAngleCmd->SetRange("Size>0.,Size<(2*pi)");
  SenAngleCmd->SetUnitCategory("Length");    
  SenAngleCmd->AvailableForStates(Idle);

  SenRotCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSenRot",this);
  SenRotCmd->SetGuidance("Set rotation angle of the sensor around z-axes");
  SenRotCmd->SetParameterName("Size",false);
  SenRotCmd->SetRange("Size>0.,Size<(2*pi)");
  SenRotCmd->SetUnitCategory("Length");    
  SenRotCmd->AvailableForStates(Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 XrayFluoDetectorMessenger::~XrayFluoDetectorMessenger()
{
  delete SamMaterCmd; delete SenMaterCmd;
  delete SamThickCmd; delete SenThickCmd;
  delete SamSizeYZCmd; delete SenSizeYZCmd;   
  delete SenDistCmd;   delete SenRotCmd;
  delete SenAngleCmd;
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == SamMaterCmd )
   { Detector->SetSampleMaterial(newValue);}
   
  if( command == SenMaterCmd )
   { Detector->SetSensorMaterial(newValue);}

  if( command == SenThickCmd )
   { Detector->SetSensorThickness(SenThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == SamThickCmd )
   { Detector->SetSampleThickness(SamThickCmd->GetNewDoubleValue(newValue));}

  if( command == SenSizeYZCmd )
   { Detector->SetSensorSizeYZ(SenSizeYZCmd->GetNewDoubleValue(newValue));}

  if( command == SamSizeYZCmd )
   { Detector->SetSampleSizeYZ(SamSizeYZCmd->GetNewDoubleValue(newValue));}
 
  if( command == SenDistCmd )
   { Detector->SetSensorDistance(SenDistCmd->GetNewDoubleValue(newValue));}
 
 if( command == SenRotCmd )
   { Detector->SetSensorAzimuth(SenRotCmd->GetNewDoubleValue(newValue));}

 if( command == SenAngleCmd )
   { Detector->SetSensorRotation(SenAngleCmd->GetNewDoubleValue(newValue));}  

 if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
