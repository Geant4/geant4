// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07DetectorMessenger.cc,v 1.2 1999-12-15 14:54:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "T07DetectorMessenger.hh"

#include "T07DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07DetectorMessenger::T07DetectorMessenger(T07DetectorConstruction * T07Det)
:T07Detector(T07Det)
{ 
  T07detDir = new G4UIdirectory("/calor/");
  T07detDir->SetGuidance("T07 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/calor/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",true);
  AbsMaterCmd->SetDefaultValue("Lead");
  AbsMaterCmd->SetCandidates("Air Aluminium Copper Galactic Lead Water");
  AbsMaterCmd->AvailableForStates(Idle);
  
  GapMaterCmd = new G4UIcmdWithAString("/calor/setGapMat",this);
  GapMaterCmd->SetGuidance("Select Material of the Gap.");
  GapMaterCmd->SetParameterName("choice",true);
  GapMaterCmd->SetDefaultValue("Air");
  GapMaterCmd->SetCandidates("Air Aluminium liquidArgon Scintillator quartz Aerogel CarbonicGas WaterSteam Galactic Beam");
  GapMaterCmd->AvailableForStates(Idle);
    
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("Size",false,false);
  AbsThickCmd->SetDefaultUnit("mm");
  AbsThickCmd->SetRange("Size>0.");
  AbsThickCmd->AvailableForStates(Idle);
  
  GapThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setGapThick",this);
  GapThickCmd->SetGuidance("Set Thickness of the Gap");
  GapThickCmd->SetParameterName("Size",false,false);
  GapThickCmd->SetDefaultUnit("mm");
  GapThickCmd->SetRange("Size>0.");
  GapThickCmd->AvailableForStates(Idle);
  
  SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setSizeYZ",this);
  SizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  SizeYZCmd->SetParameterName("Size",false,false);
  SizeYZCmd->SetDefaultUnit("mm");
  SizeYZCmd->SetRange("Size>0.");
  SizeYZCmd->AvailableForStates(Idle);
  
  NbLayersCmd = new G4UIcmdWithAnInteger("/calor/setNbOfLayers",this);
  NbLayersCmd->SetGuidance("Set number of layers.");
  NbLayersCmd->SetParameterName("NbLayers",false);
  NbLayersCmd->SetRange("NbLayers>0 && NbLayers<500");
  NbLayersCmd->AvailableForStates(Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/calor/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false,false);
  MagFieldCmd->SetDefaultUnit("tesla");
  MagFieldCmd->AvailableForStates(Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07DetectorMessenger::~T07DetectorMessenger()
{
  delete NbLayersCmd;
  delete AbsMaterCmd; delete GapMaterCmd;
  delete AbsThickCmd; delete GapThickCmd;
  delete SizeYZCmd;   delete UpdateCmd;
  delete MagFieldCmd;
  delete T07detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
   { T07Detector->SetAbsorberMaterial(newValue);}
   
  if( command == GapMaterCmd )
   { T07Detector->SetGapMaterial(newValue);}
  
  if( command == AbsThickCmd )
   { T07Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == GapThickCmd )
   { T07Detector->SetGapThickness(GapThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == SizeYZCmd )
   { T07Detector->SetCalorSizeYZ(SizeYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == NbLayersCmd )
   { T07Detector->SetNbOfLayers(NbLayersCmd->GetNewIntValue(newValue));}
  
  if( command == UpdateCmd )
   { T07Detector->UpdateGeometry(); }

  if( command == MagFieldCmd )
   { T07Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
