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
// $Id: Tst50PrimarygeneratorMessenger.cc
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "Tst50PrimaryGeneratorMessenger.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PrimaryGeneratorMessenger::Tst50PrimaryGeneratorMessenger(Tst50PrimaryGeneratorAction* Tst50Gun)
:Tst50Action(Tst50Gun)
{ 
  
  RndmDir = new G4UIcmdWithAString("/gun/randomDirection",this);
  RndmDir->SetGuidance("Shoot  incident particle with random direction.");
  RndmDir->SetGuidance("  Choice : on, off (default)");
  RndmDir->SetParameterName("choice",true);
  RndmDir->SetDefaultValue("off");
  RndmDir->SetCandidates("on off");
  RndmDir->AvailableForStates(G4State_PreInit,G4State_Idle);

  spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(G4State_PreInit,G4State_Idle);

  isoDir = new G4UIcmdWithAString("/gun/isoDir",this);
  isoDir->SetGuidance("Shoot the incident particle from an isotrofic direction.");
  isoDir->SetGuidance("  Choice : on(default), off");
  isoDir->SetParameterName("choice",true);
  isoDir->SetDefaultValue("on");
  isoDir->SetCandidates("on off");
  isoDir->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PrimaryGeneratorMessenger::~Tst50PrimaryGeneratorMessenger()
{
  
  delete  RndmDir;
  delete spectrum;
  delete isoDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst50PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmDir )
   { Tst50Action->SetRndmDirection(newValue);} 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

