// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8PhysicsListMessenger.cc,v 1.1 2000-01-07 14:50:46 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em8PhysicsListMessenger.hh"

#include "Em8PhysicsList.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em8PhysicsListMessenger::Em8PhysicsListMessenger(Em8PhysicsList * List)
:Em8List(List)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("range",true);
  cutGCmd->SetDefaultValue(1.);
  cutGCmd->SetDefaultUnit("mm");
  cutGCmd->AvailableForStates(Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/calor/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutP",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others.");
  cutPCmd->SetParameterName("range",true);
  cutPCmd->SetDefaultValue(1.);
  cutPCmd->SetDefaultUnit("mm");
  cutPCmd->AvailableForStates(Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  eCmd->SetParameterName("energy",true);
  eCmd->SetDefaultValue(10.);
  eCmd->SetDefaultUnit("keV");
  eCmd->AvailableForStates(Idle);

  rCmd = new G4UIcmdWithADoubleAndUnit("/calor/range",this);
  rCmd->SetGuidance("Display the RANGE of Electron for the current material.");
  rCmd->SetParameterName("energy",true);
  rCmd->SetDefaultValue(10.);
  rCmd->SetDefaultUnit("keV");
  rCmd->AvailableForStates(Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  setMaxStepCmd->SetGuidance("Set max. step length in the detector");
  setMaxStepCmd->SetParameterName("mxStep",true);
  setMaxStepCmd->SetDefaultUnit("mm");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em8PhysicsListMessenger::~Em8PhysicsListMessenger()
{

  delete setMaxStepCmd;

  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete rCmd;
  delete eCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void Em8PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGCmd)
    { Em8List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    { Em8List->SetElectronCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == cutPCmd)
    { Em8List->SetProtonCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == eCmd)
    { Em8List->SetCutsByEnergy(cutECmd->GetNewDoubleValue(newValue));}
  if(command == rCmd)
    { Em8List->GetRange(rCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { Em8List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

