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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PhysicsListMessenger::Tst50PhysicsListMessenger(Tst50PhysicsList * List)
:Tst50List(List)
{

  lowEnDir = new G4UIdirectory("/le/");
  lowEnDir->SetGuidance("LowEnergy commands");
 
  
  RangeDir = new G4UIcmdWithAString("/le/range_processes",this);
  RangeDir->SetGuidance("Multiple scattering and energy loss fluctuations switched.");
  RangeDir->SetGuidance("  Choice : on (default), off");
  RangeDir->SetParameterName("choice",true);
  RangeDir->SetDefaultValue("on");
  RangeDir->SetCandidates("on off");
  RangeDir->AvailableForStates(PreInit,Idle); 

 cutECmd = new G4UIcmdWithADoubleAndUnit("/le/cutE",this);
  cutECmd->SetGuidance("Set cut values by RANGE for e- e+.");
  cutECmd->SetParameterName("range",true);
  cutECmd->SetDefaultValue(1.);
  cutECmd->SetDefaultUnit("mm");
  cutECmd->AvailableForStates(Idle);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PhysicsListMessenger::~Tst50PhysicsListMessenger()
{

  

  delete  cutECmd;
  delete RangeDir;
  delete lowEnDir;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void Tst50PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  /*
  if(command == cutELowLimCmd)
    { Tst50List->SetElectronLowLimit(cutELowLimCmd->GetNewDoubleValue(newValue));}
  */
 if( command == RangeDir )
   { Tst50List->SetRangeConditions(newValue);
   G4cout<<"arrivo al PhysicsMessenger"<<G4endl;} 

if( command == cutECmd  )
 { Tst50List->SetElectronCut(cutECmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







