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
/// \file electromagnetic/TestEm10/src/Em10DetectorMessenger.cc
/// \brief Implementation of the Em10DetectorMessenger class
//
//
// $Id: Em10DetectorMessenger.cc 67268 2013-02-13 11:38:40Z ihrivnac $
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
:G4UImessenger(),
 Em10Detector(Em10Det)
{ 
  Em10detDir = new G4UIdirectory("/XTRdetector/");
  Em10detDir->SetGuidance("Em10 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/XTRdetector/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",true);
  AbsMaterCmd->SetDefaultValue("Xe");
  AbsMaterCmd->AvailableForStates(G4State_Idle);

  RadiatorMaterCmd = new G4UIcmdWithAString("/XTRdetector/setRadMat",this);
  RadiatorMaterCmd->SetGuidance("Select Material of the XTR radiator.");
  RadiatorMaterCmd->SetParameterName("choice",true);
  RadiatorMaterCmd->SetDefaultValue("CH2");
  RadiatorMaterCmd->AvailableForStates(G4State_Idle);

  DetectorSetUpCmd = new G4UIcmdWithAString("/XTRdetector/setup",this);
  DetectorSetUpCmd->SetGuidance("Select setup for comparison with experiment");
  DetectorSetUpCmd->SetParameterName("choice",true);
  DetectorSetUpCmd->SetDefaultValue("simpleALICE");
  DetectorSetUpCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ModelCmd = new G4UIcmdWithAnInteger("/XTRdetector/setModel",this);
  ModelCmd->SetGuidance("Select Model for XTR");
  ModelCmd->SetParameterName("choice",true);
  ModelCmd->SetDefaultValue(0);
  // ModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  ModelCmd->AvailableForStates(G4State_Idle);

  FoilNumCmd = new G4UIcmdWithAnInteger("/XTRdetector/setFoilNum",this);
  FoilNumCmd->SetGuidance("Select foil number for XTR");
  FoilNumCmd->SetParameterName("choice",true);
  FoilNumCmd->SetDefaultValue(0);
  FoilNumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/XTRdetector/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",true);
  WorldMaterCmd->SetDefaultValue("Air");
  WorldMaterCmd->AvailableForStates(G4State_Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false,false);
  AbsThickCmd->SetDefaultUnit("mm");
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->AvailableForStates(G4State_Idle);

  RadiatorThickCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setRadThick",this);
  RadiatorThickCmd->SetGuidance("Set Thickness of XTR radiator");
  RadiatorThickCmd->SetParameterName("SizeZ",false,false);
  RadiatorThickCmd->SetDefaultUnit("mm");
  RadiatorThickCmd->SetRange("SizeZ>0.");
  RadiatorThickCmd->AvailableForStates(G4State_Idle);

  GasGapThickCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setGasGapThick",this);
  GasGapThickCmd->SetGuidance("Set Thickness of XTR gas gaps");
  GasGapThickCmd->SetParameterName("SizeZ",false,false);
  GasGapThickCmd->SetDefaultUnit("mm");
  GasGapThickCmd->SetRange("SizeZ>0.");
  GasGapThickCmd->AvailableForStates(G4State_Idle);
  
  AbsRadCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setAbsRad",this);
  AbsRadCmd->SetGuidance("Set radius of the Absorber");
  AbsRadCmd->SetParameterName("SizeR",false,false);
  AbsRadCmd->SetDefaultUnit("mm");
  AbsRadCmd->SetRange("SizeR>0.");
  AbsRadCmd->AvailableForStates(G4State_Idle);
  
  AbsZposCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setAbsZpos",this);
  AbsZposCmd->SetGuidance("Set Z pos. of the Absorber");
  AbsZposCmd->SetParameterName("Zpos",false,false);
  AbsZposCmd->SetDefaultUnit("mm");
  AbsZposCmd->AvailableForStates(G4State_Idle);
  
  WorldZCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setWorldZ",this);
  WorldZCmd->SetGuidance("Set Z size of the World");
  WorldZCmd->SetParameterName("WSizeZ",false,false);
  WorldZCmd->SetDefaultUnit("mm");
  WorldZCmd->SetRange("WSizeZ>0.");
  WorldZCmd->AvailableForStates(G4State_Idle);
  
  WorldRCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setWorldR",this);
  WorldRCmd->SetGuidance("Set R size of the World");
  WorldRCmd->SetParameterName("WSizeR",false,false);
  WorldRCmd->SetDefaultUnit("mm");
  WorldRCmd->SetRange("WSizeR>0.");
  WorldRCmd->AvailableForStates(G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/XTRdetector/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false,false);
  MagFieldCmd->SetDefaultUnit("tesla");
  MagFieldCmd->AvailableForStates(G4State_Idle); 
  /* 
  ElectronCutCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setElectronCut",this);
  ElectronCutCmd->SetGuidance("Set electron cut in mm for vertex region");
  ElectronCutCmd->SetParameterName("ElectronCut",false,false);
  ElectronCutCmd->SetDefaultUnit("mm");
  ElectronCutCmd->SetRange("ElectronCut>0.");
  ElectronCutCmd->AvailableForStates(G4State_Idle);


  PositronCutCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setPositronCut",this);
  PositronCutCmd->SetGuidance("Set positron cut in mm for vertex region");
  PositronCutCmd->SetParameterName("PositronCut",false,false);
  PositronCutCmd->SetDefaultUnit("mm");
  PositronCutCmd->SetRange("PositronCut>0.");
  PositronCutCmd->AvailableForStates(G4State_Idle);


  GammaCutCmd = new G4UIcmdWithADoubleAndUnit("/XTRdetector/setGammaCut",this);
  GammaCutCmd->SetGuidance("Set gamma cut in mm for vertex region");
  GammaCutCmd->SetParameterName("GammaCut",false,false);
  GammaCutCmd->SetDefaultUnit("mm");
  GammaCutCmd->SetRange("GammaCut>0.");
  GammaCutCmd->AvailableForStates(G4State_Idle);
  */

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
  if( command == ModelCmd && newValue == 0 )
    { 
      //  Em10Detector->SetParametrisationModel(newValue);
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
    // Em10Detector->SetParametrisationModel(ModelCmd->GetNewIntValue(newValue));
    // Em10Detector->ParametrisationModel();
  }
  if( command == FoilNumCmd )
  { 
    Em10Detector->SetFoilNumber(FoilNumCmd->GetNewIntValue(newValue));
  }
  if( command == AbsMaterCmd )
   { Em10Detector->SetAbsorberMaterial(newValue);}

  if( command == DetectorSetUpCmd )
   { Em10Detector->SetDetectorSetUp(newValue);}

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
  /*
  if( command == ElectronCutCmd )
  { 
    Em10Detector->SetElectronCut(WorldRCmd->GetNewDoubleValue(newValue));
  }
  if( command == PositronCutCmd )
  { 
    Em10Detector->SetPositronCut(WorldRCmd->GetNewDoubleValue(newValue));
  }
  if( command == GammaCutCmd )
  { 
    Em10Detector->SetGammaCut(WorldRCmd->GetNewDoubleValue(newValue));
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
