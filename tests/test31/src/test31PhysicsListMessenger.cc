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

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/cutPositron",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others");
  cutPCmd->SetParameterName("cutHadron",false);
  cutPCmd->SetRange("cutHadron>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(PreInit,Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/cutAll",this);
  eCmd->SetGuidance("Set cut values by ENERGY for charged particles.");
  eCmd->SetParameterName("cutForAll",false);
  eCmd->SetRange("cutForAll>0.");
  eCmd->SetUnitCategory("Length");   
  eCmd->AvailableForStates(PreInit,Idle);


  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/test31/physics/MaxStep",this);
  setMaxStepCmd->SetGuidance("Set max charged particle step length");
  setMaxStepCmd->SetParameterName("MaxStep",false);
  setMaxStepCmd->SetRange("MaxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");
  setMaxStepCmd->AvailableForStates(PreInit,Idle);

  EMPhysicsCmd = new G4UIcmdWithAString("/test31/physics/addPhysics",this);
  EMPhysicsCmd->SetGuidance("Add modular PhysicsList");
  EMPhysicsCmd->SetParameterName("EMList",false);
  EMPhysicsCmd->AvailableForStates(PreInit,Idle);

  tCmd = new G4UIcmdWithAString("/test31/physics/table",this);
  tCmd->SetGuidance("Define hadron stopping table");
  tCmd->SetParameterName("h_table",false);
  tCmd->AvailableForStates(PreInit,Idle);

  HadPhysicsCmd = new G4UIcmdWithAString("/test31/physics/nuclStopping",this);
  HadPhysicsCmd->SetGuidance("Set on/off for nuclear stopping");
  HadPhysicsCmd->SetParameterName("nucStopping",false);
  HadPhysicsCmd->AvailableForStates(PreInit,Idle);

  decayCmd = new G4UIcmdWithAString("/test31/physics/Barkas",this);
  decayCmd->SetGuidance("Set on/off for Barkas effect");
  decayCmd->SetParameterName("barkas",false);
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
  delete setMaxStepCmd;
  delete EMPhysicsCmd;
  delete tCmd;
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
    { test31List->SetCutForGamma(cutGCmd->GetNewDoubleValue(newValue));}
  if(com == cutECmd)
    { test31List->SetCutForElectron(cutECmd->GetNewDoubleValue(newValue));}
  if(com == cutPCmd)
    { test31List->SetCutForPositron(cutPCmd->GetNewDoubleValue(newValue));}
  if(com == eCmd)
    { test31List->SetCutForAll(eCmd->GetNewDoubleValue(newValue));}
  if(com == setMaxStepCmd)
    { test31List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
  if(com == EMPhysicsCmd)
    { test31List->AddPhysicsList(newValue);}
  if(com == tCmd)
    { test31List->SetStoppingTable(newValue);}
  if(com == HadPhysicsCmd)
    { test31List->SetNuclearStopping(newValue);}
  if(com == decayCmd)
    { test31List->SetBarkas(newValue);}
  if( com == verbCmd )
   { test31List->SetVerbose(verbCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

