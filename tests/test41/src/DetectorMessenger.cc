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
#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");
      
  Mater1Cmd = new G4UIcmdWithAString("/testem/det/setWorldMat",this);
  Mater1Cmd->SetGuidance("Select material of world 1.");
  Mater1Cmd->SetParameterName("choice1",false);
  Mater1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  Mater2Cmd = new G4UIcmdWithAString("/testem/det/setAbsorberMat",this);
  Mater2Cmd->SetGuidance("Select material of absorber.");
  Mater2Cmd->SetParameterName("choice2",false);
  Mater2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  Mater3Cmd = new G4UIcmdWithAString("/testem/setFileName",this);
  Mater3Cmd->SetGuidance("Set output file name.");
  Mater3Cmd->SetParameterName("choice3",false);
  Mater3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeX1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/beamMomentum",this);
  SizeX1Cmd->SetGuidance("Set beam momentum");
  SizeX1Cmd->SetParameterName("SizeX1",false);
  SizeX1Cmd->SetRange("SizeX1>0.");
  SizeX1Cmd->SetUnitCategory("Energy");
  SizeX1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeX2Cmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/beamSigma",this);
  SizeX2Cmd->SetGuidance("Set beam energy spread");
  SizeX2Cmd->SetParameterName("SizeX2",false);
  SizeX2Cmd->SetRange("SizeX2>0.");
  SizeX2Cmd->SetUnitCategory("Energy");
  SizeX2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeX3Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWidth",this);
  SizeX3Cmd->SetGuidance("Set absorber width in X0");
  SizeX3Cmd->SetParameterName("SizeX3",false);
  SizeX3Cmd->SetRange("SizeX3>0.");
  SizeX3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
          
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete Mater1Cmd;
  delete Mater2Cmd;
  delete Mater3Cmd;
  delete SizeX1Cmd;
  delete SizeX2Cmd;
  delete SizeX3Cmd;

  delete UpdateCmd;
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == Mater1Cmd )
   { Detector->SetWorldMaterial(newValue);}
   
  if( command == Mater2Cmd )
   { Detector->SetAbsorberMaterial(newValue);}
   
  if( command == Mater3Cmd ) {
    G4String name = newValue;
    if(name == "HISTFILENAME") {
      char* path = getenv(name);
      if (path) name = G4String(path);
      else {
        G4cout << "### PhysicsListMessenger WARNING: "
               << " environment variable HISTFILENAME is not defined"
               << G4endl;
        return;
      }
    }    
    Detector->SetFileName(name);
  }
   
  if( command == SizeX1Cmd )
   { Detector->SetBeamMomentum(SizeX1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == SizeX2Cmd )
   { Detector->SetBeamMomentumSpread(SizeX2Cmd->GetNewDoubleValue(newValue));}
   
  if( command == SizeX3Cmd )
   { Detector->SetAbsorberWidth(SizeX3Cmd->GetNewDoubleValue(newValue));}
            
  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
