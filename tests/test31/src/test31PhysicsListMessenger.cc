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
// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31PhysicsListMessenger -------
//              
//  Modified: 08.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31PhysicsListMessenger.hh"

#include "test31PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31PhysicsListMessenger::test31PhysicsListMessenger(test31PhysicsList* list)
 :test31List(list)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/cutGamma",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma");
  cutGCmd->SetParameterName("cutGamma",false);
  cutGCmd->SetRange("cutGamma>0.");
  cutGCmd->SetUnitCategory("Length");
  cutGCmd->AvailableForStates(PreInit,Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/cutElectron",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- & e+");
  cutECmd->SetParameterName("cutElectron",false);
  cutECmd->SetRange("cutElectron>0.");
  cutECmd->SetUnitCategory("Length");  
  cutECmd->AvailableForStates(PreInit,Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/cutHadron",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others");
  cutPCmd->SetParameterName("cutHadron",false);
  cutPCmd->SetRange("cutHadron>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(PreInit,Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/cutElectronEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  eCmd->SetParameterName("cutElectronEnergy",false);
  eCmd->SetRange("cutElectronEnergy>0.");
  eCmd->SetUnitCategory("Energy");   
  eCmd->AvailableForStates(PreInit,Idle);

  lowLimCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/LowLimit",this);
  lowLimCmd->SetGuidance("Set low enery limit for charged particles");
  lowLimCmd->SetParameterName("LowLimit",false);
  lowLimCmd->SetRange("LowLimit>0.");
  lowLimCmd->SetUnitCategory("Energy");   
  lowLimCmd->AvailableForStates(PreInit,Idle);

  highLimCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/HighLimit",this);
  highLimCmd->SetGuidance("Set high energy limit for charged particles");
  highLimCmd->SetParameterName("HighLimit",false);
  highLimCmd->SetRange("HighLimit>0.");
  highLimCmd->SetUnitCategory("Energy");   
  highLimCmd->AvailableForStates(PreInit,Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/MaxStep",this);
  setMaxStepCmd->SetGuidance("Set max charged particle step length");
  setMaxStepCmd->SetParameterName("MaxStep",false);
  setMaxStepCmd->SetRange("MaxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");
  setMaxStepCmd->AvailableForStates(PreInit,Idle);

  EMPhysicsCmd = new G4UIcmdWithAString("/test31/physics/EMList",this);
  EMPhysicsCmd->SetGuidance("Set the name of the EMPhysicsList");
  EMPhysicsCmd->SetParameterName("EMList",false);
  EMPhysicsCmd->AvailableForStates(PreInit,Idle);

  HadPhysicsCmd = new G4UIcmdWithAString("/test31/physics/HadronList",this);
  HadPhysicsCmd->SetGuidance("Set the name of the HadPhysicsList");
  HadPhysicsCmd->SetParameterName("HadronList",false);
  HadPhysicsCmd->AvailableForStates(PreInit,Idle);

  decayCmd = new G4UIcmdWithAString("/test31/physics/decay",this);
  decayCmd->SetGuidance("Set the name of the decayList");
  decayCmd->SetParameterName("decay",false);
  decayCmd->AvailableForStates(PreInit,Idle);

  verbCmd = new G4UIcmdWithAnInteger("/test31/physics/verbose",this);
  verbCmd->SetGuidance("Set verbose for test31");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(PreInit,Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31PhysicsListMessenger::~test31PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete eCmd;
  delete lowLimCmd;
  delete highLimCmd;
  delete setMaxStepCmd;
  delete EMPhysicsCmd;
  delete HadPhysicsCmd;
  delete decayCmd;
  delete verbCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void test31PhysicsListMessenger::SetNewValue(G4UIcommand* com, G4String newValue)
{
  if(test31List->GetVerbose() > 1) {
    G4cout << "test31PhysicsListMessenger: new value = " << newValue << G4endl;
  }

  if(com == cutGCmd)
    { test31List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(com == cutECmd)
    { test31List->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}
  if(com == cutPCmd)
    { test31List->SetProtonCut(cutPCmd->GetNewDoubleValue(newValue));}
  if(com == eCmd)
    { test31List->SetElectronCutByEnergy(eCmd->GetNewDoubleValue(newValue));}
  if(com == lowLimCmd)
    { test31List->SetLowEnergyLimit(lowLimCmd->GetNewDoubleValue(newValue));}
  if(com == highLimCmd)
    { test31List->SetHighEnergyLimit(highLimCmd->GetNewDoubleValue(newValue));}
  if(com == setMaxStepCmd)
    { test31List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
  if(com == EMPhysicsCmd)
    { test31List->SetEMPhysicsList(newValue);}
  if(com == HadPhysicsCmd)
    { test31List->SetHadronPhysicsList(newValue);}
  if(com == decayCmd)
    { test31List->SetDecay(newValue);}
  if( com == verbCmd )
   { test31List->SetVerbose(verbCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

