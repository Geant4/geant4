// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5PhysicsListMessenger.cc,v 1.1 1999-10-12 12:23:36 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em5PhysicsListMessenger.hh"

#include "Em5PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5PhysicsListMessenger::Em5PhysicsListMessenger(Em5PhysicsList * List)
:Em5List(List)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("range",false);
  cutGCmd->SetRange("range>0.");
  cutGCmd->SetUnitCategory("Length");
  cutGCmd->AvailableForStates(Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/calor/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  cutECmd->SetParameterName("range",false);
  cutECmd->SetRange("range>0.");
  cutECmd->SetUnitCategory("Length");  
  cutECmd->AvailableForStates(Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutP",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others.");
  cutPCmd->SetParameterName("range",false);
  cutPCmd->SetRange("range>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  eCmd->SetParameterName("energy",false);
  eCmd->SetRange("energy>0.");
  eCmd->SetUnitCategory("Energy");   
  eCmd->AvailableForStates(Idle);

  rCmd = new G4UIcmdWithADoubleAndUnit("/calor/range",this);
  rCmd->SetGuidance("Display the RANGE of Electron for the current material.");
  rCmd->SetParameterName("energy",false);
  rCmd->SetRange("energy>0.");
  rCmd->SetUnitCategory("Energy");     
  rCmd->AvailableForStates(Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  setMaxStepCmd->SetGuidance("Set max. step length in the detector");
  setMaxStepCmd->SetParameterName("mxStep",false);
  setMaxStepCmd->SetRange("mxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5PhysicsListMessenger::~Em5PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete rCmd;
  delete eCmd;
  delete setMaxStepCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void Em5PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGCmd)
    { Em5List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    { Em5List->SetElectronCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == cutPCmd)
    { Em5List->SetProtonCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == eCmd)
    { Em5List->SetCutsByEnergy(cutECmd->GetNewDoubleValue(newValue));}
  if(command == rCmd)
    { Em5List->GetRange(rCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { Em5List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

