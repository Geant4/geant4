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
// $Id: Em5PhysicsListMessenger.cc,v 1.6 2002-03-08 17:30:19 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em5PhysicsListMessenger.hh"

#include "Em5PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5PhysicsListMessenger::Em5PhysicsListMessenger(Em5PhysicsList * List)
:Em5List(List)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setGCut",this);
  cutGCmd->SetGuidance("Set gamma cut.");
  cutGCmd->SetParameterName("Gcut",false);
  cutGCmd->SetRange("Gcut>0.");
  cutGCmd->SetUnitCategory("Length");
  cutGCmd->AvailableForStates(PreInit,Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setECut",this);
  cutECmd->SetGuidance("Set electron cut.");
  cutECmd->SetParameterName("Ecut",false);
  cutECmd->SetRange("Ecut>0.");
  cutECmd->SetUnitCategory("Length");  
  cutECmd->AvailableForStates(PreInit,Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/setPCut",this);
  cutPCmd->SetGuidance("Set proton cut.");
  cutPCmd->SetParameterName("Pcut",false);
  cutPCmd->SetRange("Pcut>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(PreInit,Idle);

  rCmd = new G4UIcmdWithADoubleAndUnit("/run/particle/getRange",this);
  rCmd->SetGuidance("get the electron cut for the current material.");
  rCmd->SetParameterName("energy",false);
  rCmd->SetRange("energy>0.");
  rCmd->SetUnitCategory("Energy");     
  rCmd->AvailableForStates(PreInit,Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  setMaxStepCmd->SetGuidance("Set max. step length in the detector");
  setMaxStepCmd->SetParameterName("mxStep",false);
  setMaxStepCmd->SetRange("mxStep>0.");
  setMaxStepCmd->SetUnitCategory("Length");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em5PhysicsListMessenger::~Em5PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete rCmd;
  delete setMaxStepCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
void Em5PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if(command == cutGCmd)
    { Em5List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    { Em5List->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}
  if(command == cutPCmd)
    { Em5List->SetProtonCut(cutPCmd->GetNewDoubleValue(newValue));}
  if(command == rCmd)
    { Em5List->GetRange(rCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { Em5List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

