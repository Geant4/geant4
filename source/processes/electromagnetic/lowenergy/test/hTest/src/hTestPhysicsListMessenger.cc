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
//      GEANT4 hTest
//
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
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsListMessenger::hTestPhysicsListMessenger(hTestPhysicsList* list)
 :hTestList(list)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/cutGamma",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma");
  cutGCmd->SetParameterName("cutGamma",false);
  cutGCmd->SetRange("cutGamma>0.");
  cutGCmd->SetUnitCategory("Length");
  cutGCmd->AvailableForStates(PreInit,Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/cutElectron",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- & e+");
  cutECmd->SetParameterName("cutElectron",false);
  cutECmd->SetRange("cutElectron>0.");
  cutECmd->SetUnitCategory("Length");  
  cutECmd->AvailableForStates(PreInit,Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/cutHadron",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others");
  cutPCmd->SetParameterName("cutHadron",false);
  cutPCmd->SetRange("cutHadron>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(PreInit,Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/cutElectronEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  eCmd->SetParameterName("cutElectronEnergy",false);
  eCmd->SetRange("cutElectronEnergy>0.");
  eCmd->SetUnitCategory("Energy");   
  eCmd->AvailableForStates(PreInit,Idle);

  lowLimCmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/LowLimit",this);
  lowLimCmd->SetGuidance("Set low enery limit for charged particles");
  lowLimCmd->SetParameterName("LowLimit",false);
  lowLimCmd->SetRange("LowLimit>0.");
  lowLimCmd->SetUnitCategory("Energy");   
  lowLimCmd->AvailableForStates(PreInit,Idle);

  highLimCmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/HighLimit",this);
  highLimCmd->SetGuidance("Set high energy limit for charged particles");
  highLimCmd->SetParameterName("HighLimit",false);
  highLimCmd->SetRange("HighLimit>0.");
  highLimCmd->SetUnitCategory("Energy");   
  highLimCmd->AvailableForStates(PreInit,Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/hTest/physics/MaxStep",this);
  setMaxStepCmd->SetGuidance("Set max charged particle step length");
  setMaxStepCmd->SetParameterName("MaxStep",false);
  setMaxStepCmd->SetRange("MaxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");
  setMaxStepCmd->AvailableForStates(PreInit,Idle);

  EMPhysicsCmd = new G4UIcmdWithAString("/hTest/physics/EMList",this);
  EMPhysicsCmd->SetGuidance("Set the name of the EMPhysicsList");
  EMPhysicsCmd->SetParameterName("EMList",false);
  EMPhysicsCmd->AvailableForStates(PreInit,Idle);

  HadPhysicsCmd = new G4UIcmdWithAString("/hTest/physics/HadronList",this);
  HadPhysicsCmd->SetGuidance("Set the name of the HadPhysicsList");
  HadPhysicsCmd->SetParameterName("HadronList",false);
  HadPhysicsCmd->AvailableForStates(PreInit,Idle);

  verbCmd = new G4UIcmdWithAnInteger("/hTest/physics/verbose",this);
  verbCmd->SetGuidance("Set verbose for hTest");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(PreInit,Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsListMessenger::~hTestPhysicsListMessenger()
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
  delete verbCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void hTestPhysicsListMessenger::SetNewValue(G4UIcommand* com, G4String newValue)
{
  if(hTestList->GetVerbose() > 1) {
    G4cout << "hTestPhysicsListMessenger: new value = " << newValue << G4endl;
  }

  if(com == cutGCmd)
    { hTestList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(com == cutECmd)
    { hTestList->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}
  if(com == cutPCmd)
    { hTestList->SetProtonCut(cutPCmd->GetNewDoubleValue(newValue));}
  if(com == eCmd)
    { hTestList->SetElectronCutByEnergy(eCmd->GetNewDoubleValue(newValue));}
  if(com == lowLimCmd)
    { hTestList->SetLowEnergyLimit(lowLimCmd->GetNewDoubleValue(newValue));}
  if(com == highLimCmd)
    { hTestList->SetHighEnergyLimit(highLimCmd->GetNewDoubleValue(newValue));}
  if(com == setMaxStepCmd)
    { hTestList->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
  if(com == EMPhysicsCmd)
    { hTestList->SetEMPhysicsList(newValue);}
  if(com == HadPhysicsCmd)
    { hTestList->SetHadronPhysicsList(newValue);}
  if( com == verbCmd )
   { hTestList->SetVerbose(verbCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

