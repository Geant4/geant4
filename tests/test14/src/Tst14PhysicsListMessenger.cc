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
// $Id: Tst14PhysicsListMessenger.cc,v 1.10 2003-02-23 17:22:15 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Unknown (contact: Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 22 Feb 2003 MGP          Added command for building modular PhysicsList
//                          + cleaned up
//
// -------------------------------------------------------------------

#include "Tst14PhysicsListMessenger.hh"
#include "Tst14PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

Tst14PhysicsListMessenger::Tst14PhysicsListMessenger(Tst14PhysicsList * physList)
:physicsList(physList)
{
  // MGP ---- ToDo: who is responsible for deleting lowEnDir?
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
  cutSecPhotCmd->SetGuidance("Set production threshold for secondary Gamma");
  cutSecPhotCmd->SetParameterName("energy",true);
  cutSecPhotCmd->SetDefaultValue(5e-5);
  cutSecPhotCmd->SetDefaultUnit("MeV");
  cutSecPhotCmd->AvailableForStates(G4State_Idle);

  cutSecElecCmd = new G4UIcmdWithADoubleAndUnit("/lowenergy/seceleccut",this);
  cutSecElecCmd->SetGuidance("Set production threshold for secondary e-");
  cutSecElecCmd->SetParameterName("energy",true);
  cutSecElecCmd->SetDefaultValue(5e-5);
  cutSecElecCmd->SetDefaultUnit("MeV");
  cutSecElecCmd->AvailableForStates(G4State_Idle);

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

  // Auger activation
  augerCmd = new G4UIcmdWithABool("/lowenergy/auger",this);
  augerCmd->SetGuidance("Set flag Auger electrons production.");
  augerCmd->SetParameterName("Auger",true);
  augerCmd->SetDefaultValue(false);
  augerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // Building modular PhysicsList

  physicsListCmd = new G4UIcmdWithAString("/lowenergy/addPhysics",this);  
  physicsListCmd->SetGuidance("Add chunks of PhysicsList.");
  physicsListCmd->SetParameterName("physList",false);
  physicsListCmd->AvailableForStates(G4State_PreInit);  

}

Tst14PhysicsListMessenger::~Tst14PhysicsListMessenger()
{

  delete cutGLowLimCmd;
  delete cutELowLimCmd;
  delete cutGELowLimCmd;
  delete cutSecElecCmd;
  delete cutSecPhotCmd;
  delete cutGCmd;
  delete cutECmd;
  delete augerCmd;
  delete physicsListCmd;
  delete lowEnDir;
}

void Tst14PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == cutGLowLimCmd)
    { physicsList->SetGammaLowLimit(cutGLowLimCmd->GetNewDoubleValue(newValue)); }

  if (command == cutELowLimCmd)
    { physicsList->SetElectronLowLimit(cutELowLimCmd->GetNewDoubleValue(newValue)); }

  if (command == cutGELowLimCmd)
    { physicsList->SetGELowLimit(cutGELowLimCmd->GetNewDoubleValue(newValue)); }

  if (command == cutSecPhotCmd)
    { physicsList->SetLowEnSecPhotCut(cutSecPhotCmd->GetNewDoubleValue(newValue)); }

  if (command == cutSecElecCmd)
    { physicsList->SetLowEnSecElecCut(cutSecElecCmd->GetNewDoubleValue(newValue)); }

  if (command == cutGCmd)
    { physicsList->SetGammaCut(cutGCmd->GetNewDoubleValue(newValue)); }

  if (command == cutECmd)
    { physicsList->SetElectronCut(cutECmd->GetNewDoubleValue(newValue)); }

  if (command == augerCmd)
    { physicsList->ActivateAuger(augerCmd->GetNewBoolValue(newValue)); }

  if (command == physicsListCmd)
   { physicsList->AddPhysicsList(newValue); }

}






