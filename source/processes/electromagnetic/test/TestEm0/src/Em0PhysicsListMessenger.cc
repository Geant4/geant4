// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0PhysicsListMessenger.cc,v 1.1 1999-01-08 16:32:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em0PhysicsListMessenger.hh"

#include "Em0PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0PhysicsListMessenger::Em0PhysicsListMessenger(Em0PhysicsList* EvAct)
:physList(EvAct)
{ 
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("range",true);
  cutGCmd->SetDefaultValue(1.);
  cutGCmd->SetDefaultUnit("mm");
  cutGCmd->AvailableForStates(Idle);
  cutECmd = new G4UIcmdWithADoubleAndUnit("/calor/cutC",this);
  cutECmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  cutECmd->SetParameterName("energy",true);
  cutECmd->SetDefaultValue(10.);
  cutECmd->SetDefaultUnit("keV");
  cutECmd->AvailableForStates(Idle);
  rCmd = new G4UIcmdWithADoubleAndUnit("/calor/range",this);
  rCmd->SetGuidance("Display the RANGE of Electron for the current material.");
  rCmd->SetParameterName("energy",true);
  rCmd->SetDefaultValue(10.);
  rCmd->SetDefaultUnit("keV");
  rCmd->AvailableForStates(Idle);
  eCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutE",this);
  eCmd->SetGuidance("Set cut values by RANGE for electron.");
  eCmd->SetParameterName("range",true);
  eCmd->SetDefaultValue(1.);
  eCmd->SetDefaultUnit("mm");
  eCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0PhysicsListMessenger::~Em0PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete rCmd;
  delete eCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0PhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == cutGCmd)
    {physList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    {physList->SetCutsByEnergy(cutECmd->GetNewDoubleValue(newValue));}
  if(command == rCmd)
    {physList->GetRange(rCmd->GetNewDoubleValue(newValue));}
  if(command == eCmd)
    {physList->SetECut(eCmd->GetNewDoubleValue(newValue));}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
