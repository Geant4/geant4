//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: XrayFluoPhysicsListMessenger.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoPhysicsListMessenger.hh"
#include "XrayFluoPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsListMessenger::XrayFluoPhysicsListMessenger(XrayFluoPhysicsList * List)
:XrayFluoList(List)
{

  lowEnDir = new G4UIdirectory("/lowenergy/");
  lowEnDir->SetGuidance("LowEnergy commands");

  cutGLowLimCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/lowlimG",this);
  cutGLowLimCmd->SetGuidance("Set ENERGY low limit for Gamma.");
  cutGLowLimCmd->SetParameterName("energy",true);
  cutGLowLimCmd->SetDefaultValue(1e-3);
  cutGLowLimCmd->SetDefaultUnit("MeV");
  cutGLowLimCmd->AvailableForStates(G4State_Idle);

  cutELowLimCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/lowlimE",this);
  cutELowLimCmd->SetGuidance("Set ENERGY low limit for e-.");
  cutELowLimCmd->SetParameterName("energy",true);
  cutELowLimCmd->SetDefaultValue(1e-3);
  cutELowLimCmd->SetDefaultUnit("MeV");
  cutELowLimCmd->AvailableForStates(G4State_Idle);

  cutGELowLimCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/lowlimGE",this);
  cutGELowLimCmd->SetGuidance("Set ENERGY low limit for e- and Gamma.");
  cutGELowLimCmd->SetParameterName("energy",true);
  cutGELowLimCmd->SetDefaultValue(1e-3);
  cutGELowLimCmd->SetDefaultUnit("MeV");
  cutGELowLimCmd->AvailableForStates(G4State_Idle);

  cutSecPhotCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/secphotcut",this);
  cutSecPhotCmd->SetGuidance("Set production threshold for secondary Gamma.");
  cutSecPhotCmd->SetParameterName("energy",true);
  cutSecPhotCmd->SetDefaultValue(5e-5);
  cutSecPhotCmd->SetDefaultUnit("MeV");
  cutSecPhotCmd->AvailableForStates(G4State_Idle);
/*
  cutSecElecCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/seceleccut",this);
  cutSecElecCmd->SetGuidance("Set production threshold for secondary e-");
  cutSecElecCmd->SetParameterName("energy",true);
  cutSecElecCmd->SetDefaultValue(5e-5);
  cutSecElecCmd->SetDefaultUnit("MeV");
  cutSecElecCmd->AvailableForStates(G4State_Idle);
*/
  cutGCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/cutG",this);
  cutGCmd->SetGuidance("Set cut values by RANGE for Gamma.");
  cutGCmd->SetParameterName("range",true);
  cutGCmd->SetDefaultValue(1.);
  cutGCmd->SetDefaultUnit("mm");
  cutGCmd->AvailableForStates(G4State_Idle);

  cutECmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e-.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(G4State_Idle);

  cutPCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/cutP",this);
  cutPCmd->SetGuidance("Set cut values by RANGE for proton and others.");
  cutPCmd->SetParameterName("cutP",false);
  cutPCmd->SetRange("cutP>0.");
  cutPCmd->SetUnitCategory("Length");    
  cutPCmd->AvailableForStates(G4State_Idle);
 
  eCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/cutEnergy",this);
  eCmd->SetGuidance("Set cut values by ENERGY for all particles.");
  eCmd->SetParameterName("cutenergy",false);
  eCmd->SetRange("cutenergy>0.");
  eCmd->SetUnitCategory("Energy");   
  eCmd->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsListMessenger::~XrayFluoPhysicsListMessenger()
{

  delete cutGLowLimCmd;
  delete cutELowLimCmd;
  delete cutGELowLimCmd;
  //delete cutSecElecCmd;
  delete cutSecPhotCmd;
  delete cutGCmd;
  delete cutECmd;
  delete cutPCmd;
  delete eCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void XrayFluoPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
   if(command == cutGLowLimCmd)
     { XrayFluoList->SetGammaLowLimit(cutGLowLimCmd->GetNewDoubleValue(newValue));}

   if(command == cutELowLimCmd)
     { XrayFluoList->SetElectronLowLimit(cutELowLimCmd->GetNewDoubleValue(newValue));}

   if(command == cutGELowLimCmd)
     { XrayFluoList->SetGELowLimit(cutGELowLimCmd->GetNewDoubleValue(newValue));}

  if(command == cutSecPhotCmd)
    { XrayFluoList->SetLowEnSecPhotCut(cutSecPhotCmd->GetNewDoubleValue(newValue));}

//  if(command == cutSecElecCmd)
//    { XrayFluoList->SetLowEnSecElecCut(cutSecElecCmd->GetNewDoubleValue(newValue));}

  if(command == cutGCmd)
    { XrayFluoList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue));}

  if(command == cutECmd)
    { XrayFluoList->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}

  if(command == cutPCmd)
    { XrayFluoList->SetProtonCut(cutPCmd->GetNewDoubleValue(newValue));}

  if(command == eCmd)
    { XrayFluoList->SetCutsByEnergy(eCmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






