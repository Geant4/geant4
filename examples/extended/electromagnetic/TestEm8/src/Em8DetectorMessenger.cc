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
//
// $Id: Em8DetectorMessenger.cc,v 1.7 2006/06/29 17:00:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 

#include "Em8DetectorMessenger.hh"

#include "Em8DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

///////////////////////////////////////////////////////////////////////////////////

Em8DetectorMessenger::Em8DetectorMessenger(Em8DetectorConstruction * Em8Det)
:Em8Detector(Em8Det)
{ 
  Em8detDir = new G4UIdirectory("/calor/");
  Em8detDir->SetGuidance("Em8 detector control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/calor/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",true);
  AbsMaterCmd->SetDefaultValue("Lead");
  AbsMaterCmd->AvailableForStates(G4State_Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/calor/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",true);
  WorldMaterCmd->SetDefaultValue("Air");
  WorldMaterCmd->AvailableForStates(G4State_Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false,false);
  AbsThickCmd->SetDefaultUnit("mm");
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->AvailableForStates(G4State_Idle);
  
  AbsRadCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsRad",this);
  AbsRadCmd->SetGuidance("Set radius of the Absorber");
  AbsRadCmd->SetParameterName("SizeR",false,false);
  AbsRadCmd->SetDefaultUnit("mm");
  AbsRadCmd->SetRange("SizeR>0.");
  AbsRadCmd->AvailableForStates(G4State_Idle);
  
  AbsZposCmd = new G4UIcmdWithADoubleAndUnit("/calor/setAbsZpos",this);
  AbsZposCmd->SetGuidance("Set Z pos. of the Absorber");
  AbsZposCmd->SetParameterName("Zpos",false,false);
  AbsZposCmd->SetDefaultUnit("mm");
  AbsZposCmd->AvailableForStates(G4State_Idle);
  
  WorldZCmd = new G4UIcmdWithADoubleAndUnit("/calor/setWorldZ",this);
  WorldZCmd->SetGuidance("Set Z size of the World");
  WorldZCmd->SetParameterName("WSizeZ",false,false);
  WorldZCmd->SetDefaultUnit("mm");
  WorldZCmd->SetRange("WSizeZ>0.");
  WorldZCmd->AvailableForStates(G4State_Idle);
  
  WorldRCmd = new G4UIcmdWithADoubleAndUnit("/calor/setWorldR",this);
  WorldRCmd->SetGuidance("Set R size of the World");
  WorldRCmd->SetParameterName("WSizeR",false,false);
  WorldRCmd->SetDefaultUnit("mm");
  WorldRCmd->SetRange("WSizeR>0.");
  WorldRCmd->AvailableForStates(G4State_Idle);


  ElectronCutCmd = new G4UIcmdWithADoubleAndUnit("/calor/setElectronCut",this);
  ElectronCutCmd->SetGuidance("Set electron cut in mm for vertex region");
  ElectronCutCmd->SetParameterName("ElectronCut",false,false);
  ElectronCutCmd->SetDefaultUnit("mm");
  ElectronCutCmd->SetRange("ElectronCut>0.");
  ElectronCutCmd->AvailableForStates(G4State_Idle);


  PositronCutCmd = new G4UIcmdWithADoubleAndUnit("/calor/setPositronCut",this);
  PositronCutCmd->SetGuidance("Set positron cut in mm for vertex region");
  PositronCutCmd->SetParameterName("PositronCut",false,false);
  PositronCutCmd->SetDefaultUnit("mm");
  PositronCutCmd->SetRange("PositronCut>0.");
  PositronCutCmd->AvailableForStates(G4State_Idle);


  GammaCutCmd = new G4UIcmdWithADoubleAndUnit("/calor/setGammaCut",this);
  GammaCutCmd->SetGuidance("Set gamma cut in mm for vertex region");
  GammaCutCmd->SetParameterName("GammaCut",false,false);
  GammaCutCmd->SetDefaultUnit("mm");
  GammaCutCmd->SetRange("GammaCut>0.");
  GammaCutCmd->AvailableForStates(G4State_Idle);

  
  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      

}

//////////////////////////////////////////////////////////////////////////////

Em8DetectorMessenger::~Em8DetectorMessenger()
{
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsRadCmd;  
  delete AbsZposCmd; 
  delete WorldMaterCmd;
  delete WorldZCmd;
  delete WorldRCmd;
  delete UpdateCmd;
 
  delete Em8detDir;
}

//////////////////////////////////////////////////////////////////////////////////

void Em8DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
  { 
    Em8Detector->SetAbsorberMaterial(newValue);
  } 
  if( command == WorldMaterCmd )
  { 
    Em8Detector->SetWorldMaterial(newValue);
  } 
  if( command == AbsThickCmd )
  { 
    Em8Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));
  } 
  if( command == AbsRadCmd )
  { 
    Em8Detector->SetAbsorberRadius(AbsRadCmd->GetNewDoubleValue(newValue));
  } 
  if( command == AbsZposCmd )
  { 
    Em8Detector->SetAbsorberZpos(AbsZposCmd->GetNewDoubleValue(newValue));
  }
  if( command == WorldZCmd )
  { 
    Em8Detector->SetWorldSizeZ(WorldZCmd->GetNewDoubleValue(newValue));
  } 
  if( command == WorldRCmd )
  { 
    Em8Detector->SetWorldSizeR(WorldRCmd->GetNewDoubleValue(newValue));
  }
  if( command == ElectronCutCmd )
  { 
    Em8Detector->SetElectronCut(WorldRCmd->GetNewDoubleValue(newValue));
  }
  if( command == PositronCutCmd )
  { 
    Em8Detector->SetPositronCut(WorldRCmd->GetNewDoubleValue(newValue));
  }
  if( command == GammaCutCmd )
  { 
    Em8Detector->SetGammaCut(WorldRCmd->GetNewDoubleValue(newValue));
  }
  if( command == UpdateCmd )
  { 
    Em8Detector->UpdateGeometry(); 
  }


}

/////////////////////////////////////////////////////////////////////////////////
