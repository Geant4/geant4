//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "myDetectorMessenger.hh"
#include "myDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

myDetectorMessenger::myDetectorMessenger(myDetectorConstruction * Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");
      
  SamMaterCmd = new G4UIcmdWithAString("/apparate/setSamMat",this);
  SamMaterCmd->SetGuidance("Select Material of the Sample.");
  SamMaterCmd->SetParameterName("choice",false);
  SamMaterCmd->AvailableForStates(Idle);

  SiMaterCmd = new G4UIcmdWithAString("/apparate/setSiMat",this);
  SiMaterCmd->SetGuidance("Select Material of the SiSensor.");
  SiMaterCmd->SetParameterName("choice",false);
  SiMaterCmd->AvailableForStates(Idle);

  HPGeMaterCmd = new G4UIcmdWithAString("/apparate/setHPGeMat",this);
  HPGeMaterCmd->SetGuidance("Select Material of the HPGeSensor.");
  HPGeMaterCmd->SetParameterName("choice",false);
  HPGeMaterCmd->AvailableForStates(Idle);

  Dia1MaterCmd = new G4UIcmdWithAString("/apparate/setDia1Mat",this);
  Dia1MaterCmd->SetGuidance("Select Material of the 1st diaphragm.");
  Dia1MaterCmd->SetParameterName("choice",false);
  Dia1MaterCmd->AvailableForStates(Idle);

  Dia2MaterCmd = new G4UIcmdWithAString("/apparate/setDia2Mat",this);
  Dia2MaterCmd->SetGuidance("Select Material of the 2nd diaphragm.");
  Dia2MaterCmd->SetParameterName("choice",false);
  Dia2MaterCmd->AvailableForStates(Idle);

  Dia3MaterCmd = new G4UIcmdWithAString("/apparate/setDia3Mat",this);
  Dia3MaterCmd->SetGuidance("Select Material of the 3rd diaphragm.");
  Dia3MaterCmd->SetParameterName("choice",false);
  Dia3MaterCmd->AvailableForStates(Idle);


  SamThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSamThick",this);
  SamThickCmd->SetGuidance("Set Thickness of the Sample");
  SamThickCmd->SetParameterName("Size",false);
  SamThickCmd->SetRange("Size>=0.");
  SamThickCmd->SetUnitCategory("Length");
  SamThickCmd->AvailableForStates(Idle);
  
  SiThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSiThick",this);
  SiThickCmd->SetGuidance("Set Thickness of the SiSensor");
  SiThickCmd->SetParameterName("Size",false);
  SiThickCmd->SetRange("Size>=0.");
  SiThickCmd->SetUnitCategory("Length");  
  SiThickCmd->AvailableForStates(Idle);

  HPGeThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setHPGeThick",this);
  HPGeThickCmd->SetGuidance("Set Thickness of the HPGeSensor");
  HPGeThickCmd->SetParameterName("Size",false);
  HPGeThickCmd->SetRange("Size>=0.");
  HPGeThickCmd->SetUnitCategory("Length");  
  HPGeThickCmd->AvailableForStates(Idle);

  Dia1ThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia1Thick",this);
  Dia1ThickCmd->SetGuidance("Set Thickness of the 1st diaphragm");
  Dia1ThickCmd->SetParameterName("Size",false);
  Dia1ThickCmd->SetRange("Size>=0.");
  Dia1ThickCmd->SetUnitCategory("Length");  
  Dia1ThickCmd->AvailableForStates(Idle);

  Dia2ThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia2Thick",this);
  Dia2ThickCmd->SetGuidance("Set Thickness of the 2nd diaphragm");
  Dia2ThickCmd->SetParameterName("Size",false);
  Dia2ThickCmd->SetRange("Size>=0.");
  Dia2ThickCmd->SetUnitCategory("Length");  
  Dia2ThickCmd->AvailableForStates(Idle);

  Dia3ThickCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia3Thick",this);
  Dia3ThickCmd->SetGuidance("Set Thickness of the 3rd diaphragm");
  Dia3ThickCmd->SetParameterName("Size",false);
  Dia3ThickCmd->SetRange("Size>=0.");
  Dia3ThickCmd->SetUnitCategory("Length");  
  Dia3ThickCmd->AvailableForStates(Idle);
 
  SamSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSamSizeYZ",this);
  SamSizeYZCmd->SetGuidance("Set tranverse size of the sample");
  SamSizeYZCmd->SetParameterName("Size",false);
  SamSizeYZCmd->SetRange("Size>0.");
  SamSizeYZCmd->SetUnitCategory("Length");    
  SamSizeYZCmd->AvailableForStates(Idle);

  SiSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSiSizeYZ",this);
  SiSizeYZCmd->SetGuidance("Set tranverse size of the Sisensor");
  SiSizeYZCmd->SetParameterName("Size",false);
  SiSizeYZCmd->SetRange("Size>0.");
  SiSizeYZCmd->SetUnitCategory("Length");    
  SiSizeYZCmd->AvailableForStates(Idle);

  HPGeSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setHPGeSizeYZ",this);
  HPGeSizeYZCmd->SetGuidance("Set tranverse size of the HPGesensor");
  HPGeSizeYZCmd->SetParameterName("Size",false);
  HPGeSizeYZCmd->SetRange("Size>0.");
  HPGeSizeYZCmd->SetUnitCategory("Length");    
  HPGeSizeYZCmd->AvailableForStates(Idle);

  Dia1SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia1SizeYZ",this);
  Dia1SizeYZCmd->SetGuidance("Set tranverse size of the 1st diaphragm");
  Dia1SizeYZCmd->SetParameterName("Size",false);
  Dia1SizeYZCmd->SetRange("Size>0.");
  Dia1SizeYZCmd->SetUnitCategory("Length");    
  Dia1SizeYZCmd->AvailableForStates(Idle);

 Dia2SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia2SizeYZ",this);
  Dia2SizeYZCmd->SetGuidance("Set tranverse size of the 2nd diaphragm");
  Dia2SizeYZCmd->SetParameterName("Size",false);
  Dia2SizeYZCmd->SetRange("Size>0.");
  Dia2SizeYZCmd->SetUnitCategory("Length");    
  Dia2SizeYZCmd->AvailableForStates(Idle);

 Dia3SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia3SizeYZ",this);
  Dia3SizeYZCmd->SetGuidance("Set tranverse size of the 3rd diaphragm");
  Dia3SizeYZCmd->SetParameterName("Size",false);
  Dia3SizeYZCmd->SetRange("Size>0.");
  Dia3SizeYZCmd->SetUnitCategory("Length");    
  Dia3SizeYZCmd->AvailableForStates(Idle);



  SenDistCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSenDist",this);
  SenDistCmd->SetGuidance("Set distance of the detectors");
  SenDistCmd->SetParameterName("Size",false);
  SenDistCmd->SetRange("Size>0.");
  SenDistCmd->SetUnitCategory("Length");    
  SenDistCmd->AvailableForStates(Idle);

  DiaDistCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDiaDist",this);
  DiaDistCmd->SetGuidance("Set distance of the diaphragms");
  DiaDistCmd->SetParameterName("Size",false);
  DiaDistCmd->SetRange("Size>0.");
  DiaDistCmd->SetUnitCategory("Length");    
  DiaDistCmd->AvailableForStates(Idle);

  SiAngleCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSiAngle",this);
  SiAngleCmd->SetGuidance("Set angular displacement of the Sisensor");
  SiAngleCmd->SetParameterName("Size",false);
  SiAngleCmd->SetRange("Size>0.,Size<(2*pi)");
  SiAngleCmd->SetUnitCategory("Length");    
  SiAngleCmd->AvailableForStates(Idle);

  HPGeAngleCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setHPGeAngle",this);
  HPGeAngleCmd->SetGuidance("Set angular displacement of the HPGesensor");
  HPGeAngleCmd->SetParameterName("Size",false);
  HPGeAngleCmd->SetRange("Size>0.,Size<(2*pi)");
  HPGeAngleCmd->SetUnitCategory("Length");    
  HPGeAngleCmd->AvailableForStates(Idle);

  Dia1AngleCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia1Angle",this);
  Dia1AngleCmd->SetGuidance("Set angular displacement of the 1st diaphragm");
  Dia1AngleCmd->SetParameterName("Size",false);
  Dia1AngleCmd->SetRange("Size>0.,Size<(2*pi)");
  Dia1AngleCmd->SetUnitCategory("Length");    
  Dia1AngleCmd->AvailableForStates(Idle);

  Dia2AngleCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia2Angle",this);
  Dia2AngleCmd->SetGuidance("Set angular displacement of the 2nd diaphragm");
  Dia2AngleCmd->SetParameterName("Size",false);
  Dia2AngleCmd->SetRange("Size>0.,Size<(2*pi)");
  Dia2AngleCmd->SetUnitCategory("Length");    
  Dia2AngleCmd->AvailableForStates(Idle);

  Dia3AngleCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia3Angle",this);
  Dia3AngleCmd->SetGuidance("Set angular displacement of the 3rd diaphragm");
  Dia3AngleCmd->SetParameterName("Size",false);
  Dia3AngleCmd->SetRange("Size>0.,Size<(2*pi)");
  Dia3AngleCmd->SetUnitCategory("Length");    
  Dia3AngleCmd->AvailableForStates(Idle);


  SiRotCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setSiRot",this);
  SiRotCmd->SetGuidance("Set rotation angle of the Sisensor around z-axes");
  SiRotCmd->SetParameterName("Size",false);
  SiRotCmd->SetRange("Size>0.,Size<(2*pi)");
  SiRotCmd->SetUnitCategory("Length");    
  SiRotCmd->AvailableForStates(Idle);

  HPGeRotCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setHPGeRot",this);
  HPGeRotCmd->SetGuidance("Set rotation angle of the HPGesensor around z-axes");
  HPGeRotCmd->SetParameterName("Size",false);
  HPGeRotCmd->SetRange("Size>0.,Size<(2*pi)");
  HPGeRotCmd->SetUnitCategory("Length");    
  HPGeRotCmd->AvailableForStates(Idle);

  Dia1RotCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia1Rot",this);
  HPGeRotCmd->SetGuidance("Set rotation angle of the 1st diaphragm around z-axes");
  Dia1RotCmd->SetParameterName("Size",false);
  Dia1RotCmd->SetRange("Size>0.,Size<(2*pi)");
  Dia1RotCmd->SetUnitCategory("Length");    
  Dia1RotCmd->AvailableForStates(Idle);

  Dia2RotCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia2Rot",this);
  HPGeRotCmd->SetGuidance("Set rotation angle of the 2nd diaphragm around z-axes");
  Dia2RotCmd->SetParameterName("Size",false);
  Dia2RotCmd->SetRange("Size>0.,Size<(2*pi)");
  Dia2RotCmd->SetUnitCategory("Length");    
  Dia2RotCmd->AvailableForStates(Idle);

  Dia3RotCmd = new G4UIcmdWithADoubleAndUnit("/apparate/setDia3Rot",this);
  HPGeRotCmd->SetGuidance("Set rotation angle of the 3rd diaphragm around z-axes");
  Dia3RotCmd->SetParameterName("Size",false);
  Dia3RotCmd->SetRange("Size>0.,Size<(2*pi)");
  Dia3RotCmd->SetUnitCategory("Length");    
  Dia3RotCmd->AvailableForStates(Idle);
 
  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 myDetectorMessenger::~myDetectorMessenger()
{
  delete SamMaterCmd;   delete SiMaterCmd;
  delete HPGeMaterCmd;  delete Dia1MaterCmd;
  delete Dia2MaterCmd;  delete Dia3MaterCmd;

  delete SamThickCmd;   delete SiThickCmd;
  delete HPGeThickCmd;  delete Dia1ThickCmd;
  delete Dia2ThickCmd;  delete Dia3ThickCmd;

  delete SamSizeYZCmd;  delete SiSizeYZCmd;   
  delete HPGeSizeYZCmd; delete Dia1SizeYZCmd; 
  delete Dia2SizeYZCmd; delete Dia3SizeYZCmd; 

  delete SenDistCmd;    delete DiaDistCmd;

  delete SiRotCmd;      delete HPGeRotCmd;
  delete Dia1RotCmd;    delete Dia2RotCmd;
  delete Dia3RotCmd; 
 
  delete SiAngleCmd;    delete HPGeRotCmd;
  delete Dia1AngleCmd;    delete Dia2RotCmd;
  delete Dia3AngleCmd; 

  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == SamMaterCmd )
   { Detector->SetSampleMaterial(newValue);}
   
  if( command == SiMaterCmd )
   { Detector->SetSiMaterial(newValue);}

  if( command == HPGeMaterCmd )
   { Detector->SetHPGeMaterial(newValue);}

  if( command == Dia1MaterCmd )
   { Detector->SetDia1Material(newValue);}

  if( command == Dia2MaterCmd )
   { Detector->SetDia2Material(newValue);}

  if( command == Dia3MaterCmd )
   { Detector->SetDia3Material(newValue);}

  if( command == SiThickCmd )
   { Detector->SetSiThickness(SiThickCmd->GetNewDoubleValue(newValue));}
   
 if( command == HPGeThickCmd )
   { Detector->SetHPGeThickness(HPGeThickCmd->GetNewDoubleValue(newValue));}

 if( command == Dia1ThickCmd )
   { Detector->SetDia1Thickness(Dia1ThickCmd->GetNewDoubleValue(newValue));}

 if( command == Dia2ThickCmd )
   { Detector->SetDia2Thickness(Dia2ThickCmd->GetNewDoubleValue(newValue));}

 if( command == Dia3ThickCmd )
   { Detector->SetDia3Thickness(Dia3ThickCmd->GetNewDoubleValue(newValue));}

  if( command == SamThickCmd )
   { Detector->SetSampleThickness(SamThickCmd->GetNewDoubleValue(newValue));}

  if( command == SiSizeYZCmd )
   { Detector->SetSiSizeYZ(SiSizeYZCmd->GetNewDoubleValue(newValue));}

  if( command == SamSizeYZCmd )
   { Detector->SetSampleSizeYZ(SamSizeYZCmd->GetNewDoubleValue(newValue));}
 
 if( command == HPGeSizeYZCmd )
   { Detector->SetHPGeSizeYZ(SamSizeYZCmd->GetNewDoubleValue(newValue));}

   if( command == Dia1SizeYZCmd )
     { Detector->SetDia1SizeYZ(Dia1SizeYZCmd->GetNewDoubleValue(newValue));}

 if( command == Dia2SizeYZCmd )
   { Detector->SetDia2SizeYZ(Dia2SizeYZCmd->GetNewDoubleValue(newValue));}

   if( command == Dia3SizeYZCmd )
     { Detector->SetDia3SizeYZ(Dia3SizeYZCmd->GetNewDoubleValue(newValue));}


  if( command == DiaDistCmd )
   { Detector->SetDiaphragmDistance(DiaDistCmd->GetNewDoubleValue(newValue));}

  if( command == SenDistCmd )
   { Detector->SetSensorDistance(SenDistCmd->GetNewDoubleValue(newValue));}
 
 if( command == SiRotCmd )
   { Detector->SetSiAzimuth(SiRotCmd->GetNewDoubleValue(newValue));}

 if( command == HPGeRotCmd )
   { Detector->SetHPGeAzimuth(HPGeRotCmd->GetNewDoubleValue(newValue));}

 if( command == Dia1RotCmd )
   { Detector->SetDia1Azimuth(Dia1RotCmd->GetNewDoubleValue(newValue));}

 if( command == Dia2RotCmd )
   { Detector->SetDia2Azimuth(Dia2RotCmd->GetNewDoubleValue(newValue));}

 if( command == Dia3RotCmd )
   { Detector->SetDia3Azimuth(Dia3RotCmd->GetNewDoubleValue(newValue));}
 
 if( command == SiAngleCmd )
   { Detector->SetSiRotation(SiAngleCmd->GetNewDoubleValue(newValue));} 

 if( command == HPGeAngleCmd )
   { Detector->SetHPGeRotation(HPGeAngleCmd->GetNewDoubleValue(newValue));} 

 if( command == Dia1AngleCmd )
   { Detector->SetDia1Rotation(Dia1AngleCmd->GetNewDoubleValue(newValue));} 

 if( command == Dia2AngleCmd )
   { Detector->SetDia2Rotation(Dia2AngleCmd->GetNewDoubleValue(newValue));} 

 if( command == Dia3AngleCmd )
   { Detector->SetDia3Rotation(Dia3AngleCmd->GetNewDoubleValue(newValue));} 

 if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


