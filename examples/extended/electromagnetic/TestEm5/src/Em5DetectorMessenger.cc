// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5DetectorMessenger.cc,v 1.1 1999-10-12 12:23:35 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em5DetectorMessenger.hh"

#include "Em5DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5DetectorMessenger::Em5DetectorMessenger(Em5DetectorConstruction * Em5Det)
:Em5Detector(Em5Det)
{ 
  Em5detDir = new G4UIdirectory("/calor/");
  Em5detDir->SetGuidance("Em5 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/calor/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/calor/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",false);
  WorldMaterCmd->AvailableForStates(Idle);
  
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

Em5DetectorMessenger::~Em5DetectorMessenger()
{
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsSizYZCmd;  
  delete AbsXposCmd; 
  delete WorldMaterCmd;
  delete WorldXCmd;
  delete WorldYZCmd;
  delete UpdateCmd;
  delete MagFieldCmd;
  delete Em5detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
   { Em5Detector->SetAbsorberMaterial(newValue);}
   
  if( command == WorldMaterCmd )
   { Em5Detector->SetWorldMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { Em5Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == AbsSizYZCmd )
   { Em5Detector->SetAbsorberSizeYZ(AbsSizYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == AbsXposCmd )
   { Em5Detector->SetAbsorberXpos(AbsXposCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldXCmd )
   { Em5Detector->SetWorldSizeX(WorldXCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldYZCmd )
   { Em5Detector->SetWorldSizeYZ(WorldYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { Em5Detector->UpdateGeometry(); }

  if( command == MagFieldCmd )
   { Em5Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
