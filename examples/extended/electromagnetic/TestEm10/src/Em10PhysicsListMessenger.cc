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
// $Id: Em10PhysicsListMessenger.cc,v 1.5 2005-02-01 09:37:46 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em10PhysicsListMessenger.hh"

#include "Em10PhysicsList.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10PhysicsListMessenger::Em10PhysicsListMessenger(Em10PhysicsList* List)
:Em10List(List)
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

  eMinEnergyCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/eMinEnergy",this);
  eMinEnergyCmd->SetGuidance("Set cut values by energy in Photo-Comp for e-");
  eMinEnergyCmd->SetParameterName("range",true);
  eMinEnergyCmd->SetDefaultValue(1.);
  eMinEnergyCmd->SetDefaultUnit("keV");
  eMinEnergyCmd->AvailableForStates(G4State_Idle);

  gMinEnergyCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/gMinEnergy",this);
  gMinEnergyCmd->SetGuidance("Set cut values by energy in Compton for gamma");
  gMinEnergyCmd->SetParameterName("range",true);
  gMinEnergyCmd->SetDefaultValue(1.);
  gMinEnergyCmd->SetDefaultUnit("keV");
  gMinEnergyCmd->AvailableForStates(G4State_Idle);

  setMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/step/setMaxStep",this);
  setMaxStepCmd->SetGuidance("Set max. step length in the detector");
  setMaxStepCmd->SetParameterName("mxStep",true);
  setMaxStepCmd->SetDefaultUnit("mm");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10PhysicsListMessenger::~Em10PhysicsListMessenger()
{

  delete setMaxStepCmd;

  delete cutGCmd;
  delete cutECmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGCmd)
    { Em10List->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}
  if(command == cutECmd)
    { Em10List->SetElectronCut(eCmd->GetNewDoubleValue(newValue));}
  if(command == setMaxStepCmd)
    { Em10List->SetMaxStep(setMaxStepCmd->GetNewDoubleValue(newValue));}
  if(command == eMinEnergyCmd)
    { Em10List->SetMinElectronEnergy(eMinEnergyCmd->GetNewDoubleValue(newValue));}
  if(command == gMinEnergyCmd)
    { Em10List->SetMinGammaEnergy(gMinEnergyCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

