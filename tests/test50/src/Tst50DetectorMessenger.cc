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
//
// $Id: Tst50DetectorMessenger.cc,v 1.5 2003-01-16 14:11:51 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50DetectorMessenger.hh"

#include "Tst50DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50DetectorMessenger::Tst50DetectorMessenger(
                                           Tst50DetectorConstruction* Tst50Det)
:Tst50Detector(Tst50Det)
{ 
  Tst50detDir = new G4UIdirectory("/target/");
  Tst50detDir->SetGuidance("Tst50 target control.");
      
  AbsMaterCmd = new G4UIcmdWithAString("/target/setTargetMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Target.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(Idle);
  
     
  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/target/setTargetThick",this);
  AbsThickCmd->SetGuidance("Set Thickness of the target");
  AbsThickCmd->SetParameterName("Size",false);
  AbsThickCmd->SetRange("Size>=0.");
  AbsThickCmd->SetUnitCategory("Length");
  AbsThickCmd->AvailableForStates(Idle);
  
  XThickCmd = new G4UIcmdWithADoubleAndUnit("/target/setTargetX",this);
  XThickCmd->SetGuidance("Set X dimension of the target");
  XThickCmd->SetParameterName("Size",false);
  XThickCmd->SetRange("Size>=0.");
  XThickCmd->SetUnitCategory("Length");
  XThickCmd->AvailableForStates(Idle);
  
  YThickCmd = new G4UIcmdWithADoubleAndUnit("/target/setTargetY",this);
  YThickCmd->SetGuidance("Set Y dimension of the target");
  YThickCmd->SetParameterName("Size",false);
  YThickCmd->SetRange("Size>=0.");
  YThickCmd->SetUnitCategory("Length");
  YThickCmd->AvailableForStates(Idle);
 
 SetStepCmd = new G4UIcmdWithADoubleAndUnit("/target/setMaxStep",this);
 SetStepCmd->SetGuidance("Set the particle Max Step in the target ");
 SetStepCmd->SetParameterName("Size",false);
 SetStepCmd->SetRange("Size>=0.");
 SetStepCmd->SetUnitCategory("Length");
 SetStepCmd->AvailableForStates(Idle);
  
UpdateCmd = new G4UIcmdWithoutParameter("/target/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
      
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50DetectorMessenger::~Tst50DetectorMessenger()
{
 
  delete AbsMaterCmd; 
  delete AbsThickCmd; 
 delete YThickCmd; 
 delete XThickCmd; 
  delete UpdateCmd;
  delete  SetStepCmd;
    delete Tst50detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
   { Tst50Detector->SetTargetMaterial(newValue);}
   
 
  
  if( command == AbsThickCmd )
   { Tst50Detector->SetTargetThickness(AbsThickCmd
                                               ->GetNewDoubleValue(newValue));}
  
 if( command == XThickCmd )
   { Tst50Detector->SetTargetX(XThickCmd->GetNewDoubleValue(newValue));}
  
 if( command == YThickCmd )
   { Tst50Detector->SetTargetY(YThickCmd->GetNewDoubleValue(newValue));}
 

 if( command ==SetStepCmd)
   {Tst50Detector-> SetMaxStepInTarget(SetStepCmd->GetNewDoubleValue(newValue));}

 if( command == UpdateCmd )
   { Tst50Detector->UpdateGeometry(); }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
