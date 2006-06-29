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
// $Id: DetectorMessenger.cc,v 1.4 2006-06-29 16:55:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance("UI commands specific to this example.");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");
      
  AbsMaterCmd = new G4UIcmdWithAString("/testem/det/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/testem/det/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World.");
  WorldMaterCmd->SetParameterName("wchoice",false);
  WorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAbsThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  AbsThickCmd->SetParameterName("SizeZ",false);  
  AbsThickCmd->SetRange("SizeZ>0.");
  AbsThickCmd->SetUnitCategory("Length");  
  AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAbsYZ",this);
  AbsSizYZCmd->SetGuidance("Set sizeYZ of the Absorber");
  AbsSizYZCmd->SetParameterName("SizeYZ",false);
  AbsSizYZCmd->SetRange("SizeYZ>0.");
  AbsSizYZCmd->SetUnitCategory("Length");
  AbsSizYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AbsXposCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setAbsXpos",this);
  AbsXposCmd->SetGuidance("Set X pos. of the Absorber");
  AbsXposCmd->SetParameterName("Xpos",false);
  AbsXposCmd->SetUnitCategory("Length");
  AbsXposCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldXCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWorldX",this);
  WorldXCmd->SetGuidance("Set X size of the World");
  WorldXCmd->SetParameterName("WSizeX",false);
  WorldXCmd->SetRange("WSizeX>0.");
  WorldXCmd->SetUnitCategory("Length");
  WorldXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWorldYZ",this);
  WorldYZCmd->SetGuidance("Set sizeYZ of the World");
  WorldYZCmd->SetParameterName("WSizeYZ",false);
  WorldYZCmd->SetRange("WSizeYZ>0.");
  WorldYZCmd->SetUnitCategory("Length");
  WorldYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
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
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == AbsMaterCmd )
   {Detector->SetAbsorberMaterial(newValue);}
   
  if ( command == WorldMaterCmd )
   {Detector->SetWorldMaterial(newValue);}
   
  if ( command == AbsThickCmd )
  {Detector->SetAbsorberThickness(AbsThickCmd->GetNewDoubleValue(newValue));}
   
  if ( command == AbsSizYZCmd )
   {Detector->SetAbsorberSizeYZ(AbsSizYZCmd->GetNewDoubleValue(newValue));}
   
  if ( command == AbsXposCmd )
   {Detector->SetAbsorberXpos(AbsXposCmd->GetNewDoubleValue(newValue));}
   
  if ( command == WorldXCmd )
   {Detector->SetWorldSizeX(WorldXCmd->GetNewDoubleValue(newValue));}
   
  if ( command == WorldYZCmd )
   {Detector->SetWorldSizeYZ(WorldYZCmd->GetNewDoubleValue(newValue));}
   
  if  ( command == UpdateCmd )
   {Detector->UpdateGeometry(); }

  if( command == MagFieldCmd )
   {Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
