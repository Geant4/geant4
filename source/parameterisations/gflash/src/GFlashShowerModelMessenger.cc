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
// $Id: GFlashShowerModelMessenger.cc,v 1.5 2005/11/28 18:09:26 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// ------------------------------------------------------------
// GEANT 4 class implementation
//
//      ------------- GFlashShowerModelMessenger -------------
//
// Author: Joanna Weng - 9.11.2004
// ------------------------------------------------------------

#include "GFlashShowerModelMessenger.hh"
#include "GFlashShowerModel.hh"
#include "GFlashParticleBounds.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh" 
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"

#include <iomanip>                
#include <sstream>

GFlashShowerModelMessenger::
GFlashShowerModelMessenger(GFlashShowerModel * aModel)
{ 
  myParaDir = new G4UIdirectory("/GFlash/");
  myParaDir->SetGuidance("Parametrisation control.");
  myModel= aModel;
  
  FlagCmd = new G4UIcmdWithAnInteger("/GFlash/flag",this);
  FlagCmd->SetGuidance("Defines if GFlash is activated");
  FlagCmd->SetParameterName("flag",false,false);
  
  ContCmd = new G4UIcmdWithAnInteger("/GFlash/containment ",this);
  ContCmd->SetGuidance("Defines if Containment is checked");
  ContCmd->SetParameterName("flag",false,false);
  
  StepInX0Cmd = new G4UIcmdWithADouble("/GFlash/stepXo",this);
  StepInX0Cmd->SetGuidance("Defines step lenghts");
  StepInX0Cmd->SetParameterName("flag",false,false);
  
  EminCmd = new G4UIcmdWithADoubleAndUnit("/GFlash/Emin",this);
  EminCmd->SetGuidance("Set minimum kinetic energy to trigger parametrisation");
  EminCmd->SetParameterName("Emin",false,false);
  EminCmd->SetDefaultUnit("GeV");
  EminCmd->SetUnitCategory("Energy");
  EminCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  EmaxCmd = new G4UIcmdWithADoubleAndUnit("/GFlash/Emax",this);
  EmaxCmd->SetGuidance("Set maximum kinetic energy to trigger parametrisation");
  EmaxCmd->SetParameterName("Emax",false,false);
  EmaxCmd->SetDefaultUnit("GeV");
  EmaxCmd->SetUnitCategory("Energy");
  EmaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  EkillCmd = new G4UIcmdWithADoubleAndUnit("/GFlash/Ekill",this);
  EkillCmd->SetGuidance("Set maximum kinetic energy for electrons to be killed");
  EkillCmd->SetParameterName("Ekill",false,false);
  EkillCmd->SetDefaultUnit("GeV");
  EkillCmd->SetUnitCategory("Energy");
  EkillCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}


GFlashShowerModelMessenger::~GFlashShowerModelMessenger()
{
  delete ContCmd;
  delete FlagCmd;
  delete StepInX0Cmd;  
  delete EminCmd;
  delete EmaxCmd;
  delete EkillCmd;
}


void GFlashShowerModelMessenger::
SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  
  if( command == FlagCmd ) { 
    myModel->SetFlagParamType(FlagCmd->GetNewIntValue(newValues));      
    this->GetCurrentValue(command);    
  }
  if( command == ContCmd ) { 
    myModel->SetFlagParticleContainment(ContCmd->GetNewIntValue(newValues));      
    this->GetCurrentValue(command);    
  }
  if( command == StepInX0Cmd ) { 
    myModel->SetStepInX0(StepInX0Cmd->GetNewDoubleValue(newValues));      
    this->GetCurrentValue(command);    
  }
  
  else if( command == EminCmd ) {
    myModel->PBound->SetMinEneToParametrise(*G4Electron::ElectronDefinition(),
                                       EminCmd->GetNewDoubleValue(newValues));
    this->GetCurrentValue(command);  
  }
  
  else if( command == EmaxCmd ) {
    myModel->PBound->SetMaxEneToParametrise(*G4Electron::ElectronDefinition(),
                                       EmaxCmd->GetNewDoubleValue(newValues));
    this->GetCurrentValue(command);      
  }
  
  else if( command == EkillCmd ) {
    myModel->PBound->SetEneToKill(*G4Electron::ElectronDefinition(),
                                       EkillCmd->GetNewDoubleValue(newValues));
    this->GetCurrentValue(command);  
  }
  
}


G4String GFlashShowerModelMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String returnValue('\0');
  std::ostringstream os;
  
  if( command == FlagCmd ) { 
    os << "/GFlash/flag " << myModel->GetFlagParamType()  << '\0';
    returnValue = G4String(os.str());
  }
  
  else if( command == EkillCmd ) {    
    os << "/GFlash/Ekill "
       << myModel->PBound->GetEneToKill(*G4Electron::ElectronDefinition())/GeV
       << " GeV" << '\0';
    returnValue = G4String(os.str());
  }
  
  else if( command == EminCmd ) {    
    os << "/GFlash/Emin "
       << myModel->PBound->GetMinEneToParametrise(*G4Electron::ElectronDefinition())/GeV
       << " GeV" << '\0';
    returnValue = G4String(os.str());  
  }
  
  else if( command == EmaxCmd ) {
    os << "/GFlash/Emax "
       << myModel->PBound->GetMaxEneToParametrise(*G4Electron::ElectronDefinition())/GeV
       << " GeV" << '\0';
    returnValue = G4String(os.str());
  }
  
  return returnValue;
}
