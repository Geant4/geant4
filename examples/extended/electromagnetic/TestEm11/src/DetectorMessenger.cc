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
// $Id: DetectorMessenger.cc,v 1.1 2005/06/03 15:20:32 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
  testemDir->SetGuidance(" detector control.");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");
      
  MaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  MaterCmd->SetGuidance("Select material of the box.");
  MaterCmd->SetParameterName("choice",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeXCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeX",this);
  SizeXCmd->SetGuidance("Set sizeX of the absorber");
  SizeXCmd->SetParameterName("SizeX",false);
  SizeXCmd->SetRange("SizeX>0.");
  SizeXCmd->SetUnitCategory("Length");
  SizeXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeYZ",this);
  SizeYZCmd->SetGuidance("Set sizeYZ of the absorber");
  SizeYZCmd->SetParameterName("SizeYZ",false);
  SizeYZCmd->SetRange("SizeYZ>0.");
  SizeYZCmd->SetUnitCategory("Length");
  SizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  NbLayersCmd = new G4UIcmdWithAnInteger("/testem/det/setNbOfLayers",this);
  NbLayersCmd->SetGuidance("Set number of layers");
  NbLayersCmd->SetParameterName("NbLayers",false);
  NbLayersCmd->SetRange("NbLayers>0");
  NbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
          
  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);  
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
      
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete SizeXCmd;
  delete SizeYZCmd;
  delete NbLayersCmd;   
  delete MagFieldCmd;
  delete UpdateCmd;
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MaterCmd )
   { Detector->SetMaterial(newValue);}
   
  if( command == SizeXCmd )
   { Detector->SetSizeX(SizeXCmd->GetNewDoubleValue(newValue));}
   
  if( command == SizeYZCmd )
   { Detector->SetSizeYZ(SizeYZCmd->GetNewDoubleValue(newValue));}
   
  if( command == NbLayersCmd )
   { Detector->SetNbOfLayers(NbLayersCmd->GetNewIntValue(newValue));}
         
  if( command == MagFieldCmd )
   { Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}
              
  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
