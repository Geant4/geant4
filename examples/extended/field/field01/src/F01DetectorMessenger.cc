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
// $Id: F01DetectorMessenger.cc,v 1.5 2006-06-29 17:16:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F01DetectorMessenger.hh"

#include "F01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01DetectorMessenger::F01DetectorMessenger(F01DetectorConstruction * F01Det)
:F01Detector(F01Det)
{ 
  F01detDir = new G4UIdirectory("/calor/");
  F01detDir->SetGuidance("F01 detector control.");
      
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
  
  UpdateCmd = new G4UIcmdWithoutParameter("/calor/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01DetectorMessenger::~F01DetectorMessenger()
{
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
  delete AbsRadCmd;  
  delete AbsZposCmd; 
  delete WorldMaterCmd;
  delete WorldZCmd;
  delete WorldRCmd;
  delete UpdateCmd;
  delete F01detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
   { F01Detector->SetAbsorberMaterial(newValue);}
   
  if( command == WorldMaterCmd )
   { F01Detector->SetWorldMaterial(newValue);}
   
  if( command == AbsThickCmd )
   { F01Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == AbsRadCmd )
   { F01Detector->SetAbsorberRadius(AbsRadCmd->GetNewDoubleValue(newValue));}
   
  if( command == AbsZposCmd )
   { F01Detector->SetAbsorberZpos(AbsZposCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldZCmd )
   { F01Detector->SetWorldSizeZ(WorldZCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldRCmd )
   { F01Detector->SetWorldSizeR(WorldRCmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { F01Detector->UpdateGeometry(); }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
