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
// $Id: XrayFluoPhysicsListMessenger.cc
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "Tst50PhysicsListMessenger.hh"
#include "Tst50PhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PhysicsListMessenger::Tst50PhysicsListMessenger(Tst50PhysicsList * List)
:Tst50List(List)
{
  //lowEnDir = new G4UIdirectory("/le/");
  //lowEnDir->SetGuidance("LowEnergy commands");
 
  EnDir = new G4UIdirectory("/physics/");
  EnDir->SetGuidance("physics commands");

  physicsListCmd = new G4UIcmdWithAString("/physics/addPhysics",this);  
  physicsListCmd->SetGuidance("Add chunks of PhysicsList:photon-standard, photon-epdl, photon-penelope");
  physicsListCmd->SetParameterName("choice",false);
  physicsListCmd->AvailableForStates(G4State_PreInit);  

  
 cutECmd = new G4UIcmdWithADoubleAndUnit("/physics/cutE",this);
  cutECmd->SetGuidance("Set cut values.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(G4State_Idle);
  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PhysicsListMessenger::~Tst50PhysicsListMessenger()
{

  
  
  delete cutECmd;
  delete  physicsListCmd;
  delete EnDir;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void Tst50PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  

if (command == physicsListCmd)
   { Tst50List->AddPhysicsList(newValue); }

 if (command == cutECmd)
   {Tst50List->SetParticleCut(cutECmd->GetNewDoubleValue(newValue)); }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







