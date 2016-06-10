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
// $Id: G4ErrorMessenger.cc 66892 2013-01-17 10:57:59Z gunter $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorMessenger.hh"
#include "G4ErrorStepLengthLimitProcess.hh"
#include "G4ErrorMagFieldLimitProcess.hh"
#include "G4ErrorEnergyLoss.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "globals.hh"

#ifdef G4VERBOSE
#include "G4ErrorPropagatorData.hh" //for verbosity checking
#endif

//------------------------------------------------------------------------
G4ErrorMessenger::G4ErrorMessenger(G4ErrorStepLengthLimitProcess* lengthAct, G4ErrorMagFieldLimitProcess* magAct, G4ErrorEnergyLoss* elossAct)
  :StepLengthAction(lengthAct),MagFieldAction(magAct),EnergyLossAction(elossAct)
{

  myDir = new G4UIdirectory("/geant4e/");
  myDir->SetGuidance("GEANT4e control commands");
  
  myDirLimits = new G4UIdirectory("/geant4e/limits/");
  myDirLimits->SetGuidance("GEANT4e commands to limit the step");


  StepLengthLimitCmd = new G4UIcmdWithADoubleAndUnit("/geant4e/limits/stepLength",this);  
  StepLengthLimitCmd->SetGuidance("Limit the length of an step");
  StepLengthLimitCmd->SetDefaultUnit("mm");
  StepLengthLimitCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  MagFieldLimitCmd = new G4UIcmdWithADouble("/geant4e/limits/magField",this);  
  MagFieldLimitCmd->SetGuidance("Limit the length of an step");
  MagFieldLimitCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  EnergyLossCmd = new G4UIcmdWithADouble("/geant4e/limits/energyLoss",this);  
  EnergyLossCmd->SetGuidance("Limit the length of an step");
  EnergyLossCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_GeomClosed,G4State_EventProc);
}


//------------------------------------------------------------------------
G4ErrorMessenger::~G4ErrorMessenger()
{
  delete StepLengthLimitCmd;
  delete MagFieldLimitCmd;
  delete EnergyLossCmd;
  delete myDir;
  delete myDirLimits;
}


//------------------------------------------------------------------------
void G4ErrorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{   
  if( command == StepLengthLimitCmd ) { 
#ifdef G4VERBOSE
    if(G4ErrorPropagatorData::verbose() >= 3 ) { 
      G4cout << " G4ErrorMessenger::StepLengthAction SetStepLimit " << StepLengthLimitCmd->GetNewDoubleValue(newValue) << G4endl;
    }
#endif
    StepLengthAction->SetStepLimit(StepLengthLimitCmd->GetNewDoubleValue(newValue));
  } else if( command == MagFieldLimitCmd ) { 
#ifdef G4VERBOSE
    if(G4ErrorPropagatorData::verbose() >= 3 ) { 
      G4cout << " G4ErrorMessenger::MagFieldAction SetStepLimit " << MagFieldLimitCmd->GetNewDoubleValue(newValue) << G4endl;
    }
#endif
    MagFieldAction->SetStepLimit(MagFieldLimitCmd->GetNewDoubleValue(newValue));
  } else if( command == EnergyLossCmd ) { 
#ifdef G4VERBOSE
    if(G4ErrorPropagatorData::verbose() >= 3 ) { 
      G4cout << " G4ErrorMessenger::EnergyLossAction SetStepLimit " << EnergyLossCmd->GetNewDoubleValue(newValue) << G4endl;
    }
#endif
    EnergyLossAction->SetStepLimit(EnergyLossCmd->GetNewDoubleValue(newValue));
  }
}

