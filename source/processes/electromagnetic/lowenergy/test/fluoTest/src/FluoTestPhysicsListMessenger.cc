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
// $Id: FluoTestPhysicsListMessenger.cc,v 1.13 2001-10-31 12:33:45 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestPhysicsListMessenger.hh"
#include "FluoTestPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPhysicsListMessenger::FluoTestPhysicsListMessenger(FluoTestPhysicsList * List)
:FluoTestList(List)
{

  lowEnDir = new G4UIdirectory("/lowenergy/");
  lowEnDir->SetGuidance("LowEnergy commands");

  cutGLowLimCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/lowlimG",this);
  cutGLowLimCmd->SetGuidance("Set ENERGY low limit for Gamma.");
  cutGLowLimCmd->SetParameterName("energy",true);
  cutGLowLimCmd->SetDefaultValue(1e-3);
  cutGLowLimCmd->SetDefaultUnit("MeV");
  cutGLowLimCmd->AvailableForStates(Idle);

  cutELowLimCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/lowlimE",this);
  cutELowLimCmd->SetGuidance("Set ENERGY low limit for e-.");
  cutELowLimCmd->SetParameterName("energy",true);
  cutELowLimCmd->SetDefaultValue(1e-3);
  cutELowLimCmd->SetDefaultUnit("MeV");
  cutELowLimCmd->AvailableForStates(Idle);

  cutGELowLimCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/lowlimGE",this);
  cutGELowLimCmd->SetGuidance("Set ENERGY low limit for e- and Gamma.");
  cutGELowLimCmd->SetParameterName("energy",true);
  cutGELowLimCmd->SetDefaultValue(1e-3);
  cutGELowLimCmd->SetDefaultUnit("MeV");
  cutGELowLimCmd->AvailableForStates(Idle);

  cutSecPhotCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/secphotcut",this);
  cutSecPhotCmd->SetGuidance("Set production threshold for secondary Gamma.");
  cutSecPhotCmd->SetParameterName("energy",true);
  cutSecPhotCmd->SetDefaultValue(5e-5);
  cutSecPhotCmd->SetDefaultUnit("MeV");
  cutSecPhotCmd->AvailableForStates(Idle);

  cutSecElecCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/seceleccut",this);
  cutSecElecCmd->SetGuidance("Set production threshold for secondary e-");
  cutSecElecCmd->SetParameterName("energy",true);
  cutSecElecCmd->SetDefaultValue(5e-5);
  cutSecElecCmd->SetDefaultUnit("MeV");
  cutSecElecCmd->AvailableForStates(Idle);

  cutGCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("range",true);
  cutGCmd->SetDefaultValue(1.);
  cutGCmd->SetDefaultUnit("mm");
  cutGCmd->AvailableForStates(Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e-.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPhysicsListMessenger::~FluoTestPhysicsListMessenger()
{

  delete cutGLowLimCmd;
  delete cutELowLimCmd;
  delete cutGELowLimCmd;
  delete cutSecElecCmd;
  delete cutSecPhotCmd;
  delete cutGCmd;
  delete cutECmd;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void FluoTestPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGLowLimCmd)
    { FluoTestList->SetGammaLowLimit(cutGLowLimCmd->GetNewDoubleValue(newValue));}

  if(command == cutELowLimCmd)
    { FluoTestList->SetElectronLowLimit(cutELowLimCmd->GetNewDoubleValue(newValue));}

  if(command == cutGELowLimCmd)
    { FluoTestList->SetGELowLimit(cutGELowLimCmd->GetNewDoubleValue(newValue));}

  if(command == cutSecPhotCmd)
    { FluoTestList->SetLowEnSecPhotCut(cutSecPhotCmd->GetNewDoubleValue(newValue));}

  if(command == cutSecElecCmd)
    { FluoTestList->SetLowEnSecElecCut(cutSecElecCmd->GetNewDoubleValue(newValue));}

  if(command == cutGCmd)
    { FluoTestList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}

  if(command == cutECmd)
    { FluoTestList->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






