// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestPhysicsListMessenger -------
//              
//  Modified: 08.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestPhysicsListMessenger.hh"

#include "hTestPhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsListMessenger::hTestPhysicsListMessenger(hTestPhysicsList* list)
 :hTestList(list)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/hTest/cutGamma",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("cutGamma",false);
  cutGCmd->SetRange("cutGamma>0.");
  cutGCmd->SetUnitCategory("Length");
  cutGCmd->AvailableForStates(Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/hTest/cutElectron",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  cutECmd->SetParameterName("cutElectron",false);
  cutECmd->SetRange("cutElectron>0.");
  cutECmd->SetUnitCategory("Length");  
  cutECmd->AvailableForStates(Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/hTest/cutHadron",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others.");
  cutPCmd->SetParameterName("cutHadron",false);
  cutPCmd->SetRange("cutHadron>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/hTest/cutElectronEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  eCmd->SetParameterName("cutElectronEnergy",false);
  eCmd->SetRange("cutElectronEnergy>0.");
  eCmd->SetUnitCategory("Energy");   
  eCmd->AvailableForStates(Idle);

  lowLimCmd = new G4UIcmdWithADoubleAndUnit("/hTest/setLowEnergyLimit",this);
  lowLimCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  lowLimCmd->SetParameterName("setLowEnergyLimit",false);
  lowLimCmd->SetRange("setLowEnergyLimit>0.");
  lowLimCmd->SetUnitCategory("Energy");   
  lowLimCmd->AvailableForStates(Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxChargedStep",this);
  setMaxStepCmd->SetGuidance("Set max charged particle step length");
  setMaxStepCmd->SetParameterName("mxStep",false);
  setMaxStepCmd->SetRange("mxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");

  EMPhysicsCmd = new G4UIcmdWithAString("/hTest/setEMPhysics",this);
  EMPhysicsCmd->SetGuidance("Set the name of the EMPhysicsList");
  EMPhysicsCmd->SetParameterName("EMPhysics",false);
  EMPhysicsCmd->AvailableForStates(Idle);

  HadPhysicsCmd = new G4UIcmdWithAString("/hTest/setHadPhysics",this);
  HadPhysicsCmd->SetGuidance("Set the name of the HadPhysicsList");
  HadPhysicsCmd->SetParameterName("HadPhysics",false);
  HadPhysicsCmd->AvailableForStates(Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsListMessenger::~hTestPhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete eCmd;
  delete lowLimCmd;
  delete setMaxStepCmd;
  delete EMPhysicsCmd;
  delete HadPhysicsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void hTestPhysicsListMessenger::SetNewValue(G4UIcommand* com, G4String newValue)
{
  if(com == cutGCmd)
    { hTestList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(com == cutECmd)
    { hTestList->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}
  if(com == cutPCmd)
    { hTestList->SetProtonCut(cutPCmd->GetNewDoubleValue(newValue));}
  if(com == eCmd)
    { hTestList->SetElectronCutByEnergy(eCmd->GetNewDoubleValue(newValue));}
  if(com == lowLimCmd)
    { hTestList->GetLowEnergyLimit(rCmd->GetNewDoubleValue(newValue));}
  if(com == setMaxStepCmd)
    { hTestList->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
  if(com == EMPhysicsCmd)
    { hTestList->SetEMPhysicsList(newValue);}
  if(com == HadPhysicsCmd)
    { hTestList->SetHadronPhysicsList(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

