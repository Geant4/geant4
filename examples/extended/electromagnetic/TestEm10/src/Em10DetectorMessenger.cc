// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10DetectorMessenger.cc,v 1.3 2001-03-23 13:51:44 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Em10DetectorMessenger.hh"

#include "Em10DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//////////////////////////////////////////////////////////////////////////////

Em10DetectorMessenger::Em10DetectorMessenger(Em10DetectorConstruction * Em10Det)
:Em10Detector(Em10Det)
{ 
  Em10detDir = new G4UIdirectory("/XTRdetector/");
  Em10detDir->SetGuidance("Em10 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/XTRdetector/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",true);
  AbsMaterCmd->SetDefaultValue("Xe");
  AbsMaterCmd->AvailableForStates(Idle);

  RadiatorMaterCmd = new G4UIcmdWithAString("/XTRdetector/setRadMat",this);
  RadiatorMaterCmd->SetGuidance("Select Material of the XTR radiator.");
  RadiatorMaterCmd->SetParameterName("choice",true);
  RadiatorMaterCmd->SetDefaultValue("CH2");
  RadiatorMaterCmd->AvailableForStates(Idle);

  ModelCmd = new G4UIcmdWithAnInteger("/XTRdetector/setModel",this);
  ModelCmd->SetGuidance("Select Model for XTR");
  ModelCmd->SetParameterName("choice",true);
  ModelCmd->SetDefaultValue(0);
  ModelCmd->AvailableForStates(PreInit,Idle);

  FoilNumCmd = new G4UIcmdWithAnInteger("/XTRdetector/setFoilNum",this);
  FoilNumCmd->SetGuidance("Select foil number for XTR");
  FoilNumCmd->SetParameterName("choice",true);
  FoilNumCmd->SetDefaultValue(0);
  FoilNumCmd->AvailableForStates(PreInit,Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/XTRdetector/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",true);
  WorldMaterCmd->SetDefaultValue("Air");
  WorldMaterCmd->AvailableForStates(Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false,false);
  AbsThickCmd->SetDefaultUnit("mm");
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->AvailableForStates(Idle);

RadiatorThickCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setRadThick",this);
  RadiatorThickCmd->SetGuidance("Set Thickness of XTR radiator");
  RadiatorThickCmd->SetParameterName("SizeZ",false,false);
  RadiatorThickCmd->SetDefaultUnit("mm");
  RadiatorThickCmd->SetRange("SizeZ>0.");
  RadiatorThickCmd->AvailableForStates(Idle);

GasGapThickCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setGasGapThick",this);
  GasGapThickCmd->SetGuidance("Set Thickness of XTR gas gaps");
  GasGapThickCmd->SetParameterName("SizeZ",false,false);
  GasGapThickCmd->SetDefaultUnit("mm");
  GasGapThickCmd->SetRange("SizeZ>0.");
  GasGapThickCmd->AvailableForStates(Idle);
  
  AbsRadCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setAbsRad",this);
  AbsRadCmd->SetGuidance("Set radius of the Absorber");
  AbsRadCmd->SetParameterName("SizeR",false,false);
  AbsRadCmd->SetDefaultUnit("mm");
  AbsRadCmd->SetRange("SizeR>0.");
  AbsRadCmd->AvailableForStates(Idle);
  
  AbsZposCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setAbsZpos",this);
  AbsZposCmd->SetGuidance("Set Z pos. of the Absorber");
  AbsZposCmd->SetParameterName("Zpos",false,false);
  AbsZposCmd->SetDefaultUnit("mm");
  AbsZposCmd->AvailableForStates(Idle);
  
  WorldZCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setWorldZ",this);
  WorldZCmd->SetGuidance("Set Z size of the World");
  WorldZCmd->SetParameterName("WSizeZ",false,false);
  WorldZCmd->SetDefaultUnit("mm");
  WorldZCmd->SetRange("WSizeZ>0.");
  WorldZCmd->AvailableForStates(Idle);
  
  WorldRCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setWorldR",this);
  WorldRCmd->SetGuidance("Set R size of the World");
  WorldRCmd->SetParameterName("WSizeR",false,false);
  WorldRCmd->SetDefaultUnit("mm");
  WorldRCmd->SetRange("WSizeR>0.");
  WorldRCmd->AvailableForStates(Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/XTRdetector/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false,false);
  MagFieldCmd->SetDefaultUnit("tesla");
  MagFieldCmd->AvailableForStates(Idle);  

}

///////////////////////////////////////////////////////////////////////////////

Em10DetectorMessenger::~Em10DetectorMessenger()
{
  delete AbsMaterCmd; 
  delete ModelCmd; 
  delete FoilNumCmd; 
  delete AbsThickCmd; 
  delete RadiatorThickCmd; 
  delete GasGapThickCmd; 
  delete AbsRadCmd;  
  delete AbsZposCmd; 
  delete WorldMaterCmd;
  delete RadiatorMaterCmd;
  delete WorldZCmd;
  delete WorldRCmd;
  delete UpdateCmd;
  delete MagFieldCmd;
  delete Em10detDir;
}

/////////////////////////////////////////////////////////////////////////////

void Em10DetectorMessenger::SetNewValue(G4UIcommand* command,G4int newValue)
{ 
  if( command == ModelCmd )
  { 
    Em10Detector->SetParametrisationModel(newValue);
    Em10Detector->Construct();
  }
}
////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == ModelCmd )
  { 
    Em10Detector->SetParametrisationModel(ModelCmd->GetNewIntValue(newValue));
  }
  if( command == FoilNumCmd )
  { 
    Em10Detector->SetFoilNumber(FoilNumCmd->GetNewIntValue(newValue));
  }
  if( command == AbsMaterCmd )
   { Em10Detector->SetAbsorberMaterial(newValue);}

  if( command == RadiatorMaterCmd )
   { Em10Detector->SetRadiatorMaterial(newValue);}
   
  if( command == WorldMaterCmd )
   { Em10Detector->SetWorldMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { Em10Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}

  if( command == RadiatorThickCmd )
   { Em10Detector->SetRadiatorThickness(RadiatorThickCmd->
                   GetNewDoubleValue(newValue));}

  if( command == GasGapThickCmd )
   { Em10Detector->SetGasGapThickness(GasGapThickCmd->
                   GetNewDoubleValue(newValue));}
   
  if( command == AbsRadCmd )
   { Em10Detector->SetAbsorberRadius(AbsRadCmd->GetNewDoubleValue(newValue));}
   
  if( command == AbsZposCmd )
   { Em10Detector->SetAbsorberZpos(AbsZposCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldZCmd )
   { Em10Detector->SetWorldSizeZ(WorldZCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldRCmd )
   { Em10Detector->SetWorldSizeR(WorldRCmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { Em10Detector->UpdateGeometry(); }

  if( command == MagFieldCmd )
   { Em10Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
