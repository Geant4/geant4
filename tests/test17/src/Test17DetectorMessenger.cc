// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17DetectorMessenger.hh"

#include "Test17DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17DetectorMessenger::Test17DetectorMessenger(Test17DetectorConstruction * Test17Det)
:Test17Detector(Test17Det)
{ 
  Test17detDir = new G4UIdirectory("/calor/");
  Test17detDir->SetGuidance("Test17 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/calor/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/calor/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",false);
  WorldMaterCmd->AvailableForStates(Idle);

  NumOfAbsCmd = new G4UIcmdWithAnInteger("/calor/setNumOfAbs",this);
  NumOfAbsCmd->SetGuidance("Set number of absorbers");
  NumOfAbsCmd->SetParameterName("Nabs",false);
  NumOfAbsCmd->AvailableForStates(Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false);  
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->SetUnitCategory("Length");  
  AbsThickCmd->AvailableForStates(Idle);
  
  AbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsYZ",this);
  AbsSizYZCmd->SetGuidance("Set sizeYZ of the Absorber");
  AbsSizYZCmd->SetParameterName("SizeYZ",false);
  AbsSizYZCmd->SetRange("SizeYZ>0.");
  AbsSizYZCmd->SetUnitCategory("Length");
  AbsSizYZCmd->AvailableForStates(Idle);
  
  AbsXposCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsXpos",this);
  AbsXposCmd->SetGuidance("Set X pos. of the Absorber");
  AbsXposCmd->SetParameterName("Xpos",false);
  AbsXposCmd->SetUnitCategory("Length");
  AbsXposCmd->AvailableForStates(Idle);
  
  WorldXCmd = new G4UIcmdWithADoubleAndUnit("/calor/setWorldX",this);
  WorldXCmd->SetGuidance("Set X size of the World");
  WorldXCmd->SetParameterName("WSizeX",false);
  WorldXCmd->SetRange("WSizeX>0.");
  WorldXCmd->SetUnitCategory("Length");
  WorldXCmd->AvailableForStates(Idle);
  
  WorldYZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setWorldYZ",this);
  WorldYZCmd->SetGuidance("Set sizeYZ of the World");
  WorldYZCmd->SetParameterName("WSizeYZ",false);
  WorldYZCmd->SetRange("WSizeYZ>0.");
  WorldYZCmd->SetUnitCategory("Length");
  WorldYZCmd->AvailableForStates(Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(Idle);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17DetectorMessenger::~Test17DetectorMessenger()
{
  delete NumOfAbsCmd; 
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsSizYZCmd;  
  delete AbsXposCmd; 
  delete WorldMaterCmd;
  delete WorldXCmd;
  delete WorldYZCmd;
  delete UpdateCmd;
  delete MagFieldCmd;
  delete Test17detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
   { Test17Detector->SetAbsorberMaterial(newValue);}
   
  if( command == WorldMaterCmd )
   { Test17Detector->SetWorldMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { Test17Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}

  if( command == NumOfAbsCmd )
   { Test17Detector->SetNumberOfAbsorbers(NumOfAbsCmd->GetNewIntValue(newValue));}
   
  if( command == AbsSizYZCmd )
   { Test17Detector->SetAbsorberSizeYZ(AbsSizYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == AbsXposCmd )
   { Test17Detector->SetAbsorberXpos(AbsXposCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldXCmd )
   { Test17Detector->SetWorldSizeX(WorldXCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldYZCmd )
   { Test17Detector->SetWorldSizeYZ(WorldYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { Test17Detector->UpdateGeometry(); }

  if( command == MagFieldCmd )
   { Test17Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
