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
// $Id: Tst14PhysicsListMessenger.cc,v 1.12 2006-06-29 21:42:18 gunter Exp $
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

  // Select angular distribution
  
  angularDistributionCmd = new G4UIcmdWithAString("/lowenergy/angular",this);  
  angularDistributionCmd->SetGuidance("Select angular distribution: tsai or 2bn or 2bs");
  angularDistributionCmd->SetParameterName("angularDistribution",false);
  angularDistributionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

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
  delete angularDistributionCmd;
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
  
  if (command == angularDistributionCmd)
    { physicsList->SetAngularDistribution(newValue); }
}






