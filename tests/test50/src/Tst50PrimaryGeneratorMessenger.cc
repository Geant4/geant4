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
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);

  RndmVert = new G4UIcmdWithAString("/gun/randomVert",this);
  RndmVert->SetGuidance("Shoot randomly the incident particle.");
  RndmVert->SetGuidance("  Choice : on(default), off");
  RndmVert->SetParameterName("choice",true);
  RndmVert->SetDefaultValue("on");
  RndmVert->SetCandidates("on off");
  RndmVert->AvailableForStates(PreInit,Idle);

  spectrum = new G4UIcmdWithAString("/gun/spectrum",this);
  spectrum->SetGuidance("Shoot the incident particle with a certain energy spectrum.");
  spectrum->SetGuidance("  Choice : on(default), off");
  spectrum->SetParameterName("choice",true);
  spectrum->SetDefaultValue("on");
  spectrum->SetCandidates("on off");
  spectrum->AvailableForStates(PreInit,Idle);

  isoVert = new G4UIcmdWithAString("/gun/isoVert",this);
  isoVert->SetGuidance("Shoot the incident particle from an isotrofic direction.");
  isoVert->SetGuidance("  Choice : on(default), off");
  isoVert->SetParameterName("choice",true);
  isoVert->SetDefaultValue("on");
  isoVert->SetCandidates("on off");
  isoVert->AvailableForStates(PreInit,Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PrimaryGeneratorMessenger::~Tst50PrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete  RndmVert;
  delete spectrum;
  delete isoVert;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst50PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

