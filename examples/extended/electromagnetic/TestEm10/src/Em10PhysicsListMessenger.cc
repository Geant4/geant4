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
// $Id: Em10PhysicsListMessenger.cc,v 1.6 2005/11/29 14:42:22 grichine Exp $
// GEANT4 tag $Name: geant4-08-00 $
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


  ElectronCutCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/setElectronCut",this);
  ElectronCutCmd->SetGuidance("Set electron cut in mm for vertex region");
  ElectronCutCmd->SetParameterName("ElectronCut",false,false);
  ElectronCutCmd->SetDefaultUnit("mm");
  ElectronCutCmd->SetRange("ElectronCut>0.");
  ElectronCutCmd->AvailableForStates(G4State_Idle);


  PositronCutCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/setPositronCut",this);
  PositronCutCmd->SetGuidance("Set positron cut in mm for vertex region");
  PositronCutCmd->SetParameterName("PositronCut",false,false);
  PositronCutCmd->SetDefaultUnit("mm");
  PositronCutCmd->SetRange("PositronCut>0.");
  PositronCutCmd->AvailableForStates(G4State_Idle);


  GammaCutCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/setGammaCut",this);
  GammaCutCmd->SetGuidance("Set gamma cut in mm for vertex region");
  GammaCutCmd->SetParameterName("GammaCut",false,false);
  GammaCutCmd->SetDefaultUnit("mm");
  GammaCutCmd->SetRange("GammaCut>0.");
  GammaCutCmd->AvailableForStates(G4State_Idle);

  RadiatorCutCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/setRadiatorCuts",this);
  RadiatorCutCmd->SetGuidance("Set radiator cut in mm for vertex region");
  RadiatorCutCmd->SetParameterName("RadiatorCuts",false,false);
  RadiatorCutCmd->SetDefaultUnit("mm");
  RadiatorCutCmd->SetRange("RadiatorCuts > 0.");
  RadiatorCutCmd->AvailableForStates(G4State_Idle);

  DetectorCutCmd = new G4UIcmdWithADoubleAndUnit("/emphyslist/setDetectorCuts",this);
  DetectorCutCmd->SetGuidance("Set radiator cut in mm for vertex region");
  DetectorCutCmd->SetParameterName("DetectorCuts",false,false);
  DetectorCutCmd->SetDefaultUnit("mm");
  DetectorCutCmd->SetRange("DetectorCuts > 0.");
  DetectorCutCmd->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10PhysicsListMessenger::~Em10PhysicsListMessenger()
{
  delete cutGCmd;
  delete cutECmd;

  delete eMinEnergyCmd;
  delete gMinEnergyCmd;

  delete setMaxStepCmd;

  delete ElectronCutCmd;
  delete PositronCutCmd;
  delete GammaCutCmd;
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

  if( command == ElectronCutCmd )
  {
    Em10List->SetRegElectronCut(ElectronCutCmd->GetNewDoubleValue(newValue));
  }
  if( command == PositronCutCmd )
  {
    Em10List->SetRegPositronCut(PositronCutCmd->GetNewDoubleValue(newValue));
  }
  if( command == GammaCutCmd )
  {
    Em10List->SetRegGammaCut(GammaCutCmd->GetNewDoubleValue(newValue));
  }
  if( command == RadiatorCutCmd )
  {
    Em10List->SetRadiatorCuts();
  }
  if( command == DetectorCutCmd )
  {
    Em10List->SetDetectorCuts();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

