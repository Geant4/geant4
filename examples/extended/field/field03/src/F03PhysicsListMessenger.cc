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
// $Id: F03PhysicsListMessenger.cc,v 1.4 2003/08/27 17:50:29 vnivanch Exp $
// GEANT4 tag $Name: fieldex-V05-02-00 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F03PhysicsListMessenger.hh"

#include "F03PhysicsList.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F03PhysicsListMessenger::F03PhysicsListMessenger(F03PhysicsList * List)
:F03List(List)
{
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/calor/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("range",true);
  cutGCmd->SetDefaultValue(1.);
  cutGCmd->SetDefaultUnit("mm");
  cutGCmd->AvailableForStates(G4State_Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/calor/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(G4State_Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  setMaxStepCmd->SetGuidance("Set max. step length in the detector");
  setMaxStepCmd->SetParameterName("mxStep",true);
  setMaxStepCmd->SetDefaultUnit("mm");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F03PhysicsListMessenger::~F03PhysicsListMessenger()
{

  delete setMaxStepCmd;

  delete cutGCmd;
  delete cutECmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void F03PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGCmd)
    { F03List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    { F03List->SetElectronCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { F03List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

