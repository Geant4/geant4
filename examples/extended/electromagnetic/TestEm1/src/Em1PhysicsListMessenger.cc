// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1PhysicsListMessenger.cc,v 1.1 1999-10-11 13:07:47 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1PhysicsListMessenger.hh"

#include "Em1PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1PhysicsListMessenger::Em1PhysicsListMessenger(Em1PhysicsList* EvAct)
:physList(EvAct)
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1PhysicsListMessenger::~Em1PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete rCmd;
  delete eCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsListMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == cutGCmd)
    {physList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    {physList->SetElectronCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == cutPCmd)
    {physList->SetProtonCut(eCmd->GetNewDoubleValue(newValue));}        
  if(command == eCmd)
    {physList->SetCutsByEnergy(cutECmd->GetNewDoubleValue(newValue));}
  if(command == rCmd)
    {physList->GetRange(rCmd->GetNewDoubleValue(newValue));}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
