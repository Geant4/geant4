// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1DetectorMessenger.cc,v 1.2 1999-12-15 14:48:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1DetectorMessenger.hh"

#include "Em1DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1DetectorMessenger::Em1DetectorMessenger(Em1DetectorConstruction * Em1Det)
:Em1Detector(Em1Det)
{ 
  Em1detDir = new G4UIdirectory("/calor/");
  Em1detDir->SetGuidance("Em1 detector control.");
      
  MaterCmd = new G4UIcmdWithAString("/calor/setMat",this);
  MaterCmd->SetGuidance("Select material of the box.");
  MaterCmd->SetParameterName("choice",false);
  MaterCmd->AvailableForStates(Idle);
  
  SizeCmd = new G4UIcmdWithADoubleAndUnit("/calor/setSize",this);
  SizeCmd->SetGuidance("Set size of the box");
  SizeCmd->SetParameterName("Size",false);
  SizeCmd->SetRange("Size>0.");
  SizeCmd->SetUnitCategory("Length");
  SizeCmd->AvailableForStates(Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(Idle);
    
  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1DetectorMessenger::~Em1DetectorMessenger()
{
  delete MaterCmd;
  delete SizeCmd; 
  delete MagFieldCmd;
  delete UpdateCmd;
  delete Em1detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MaterCmd )
   { Em1Detector->SetMaterial(newValue);}
   
  if( command == SizeCmd )
   { Em1Detector->SetSize(SizeCmd->GetNewDoubleValue(newValue));}
   
  if( command == MagFieldCmd )
   { Em1Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
     
  if( command == UpdateCmd )
   { Em1Detector->UpdateGeometry(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
