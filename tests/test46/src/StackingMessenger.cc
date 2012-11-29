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
// $Id$
//
/////////////////////////////////////////////////////////////////////////
//
// StackingMessenger
//
// Created: 31.05.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

#include "StackingMessenger.hh"

#include "StackingAction.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingMessenger::StackingMessenger(StackingAction* stack)
:stackAction(stack)
{
  killCmd = new G4UIcmdWithABool("/testhadr/KillAllSecondaries",this);
  killCmd->SetGuidance("  Choice : true false");
  killCmd->SetParameterName("choice",true);
  killCmd->SetDefaultValue(false);

  kCmd = new G4UIcmdWithAString("/testhadr/Kill", this);
  kCmd->SetGuidance("Kill secondary particles of defined type");
  kCmd->SetParameterName("ch", true);
  kCmd->SetDefaultValue("none");
  kCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  brCmd = new G4UIcommand("/process/had/setSecBiasing",this);
  brCmd->SetGuidance("Set Russian roullette.");
  brCmd->SetGuidance("  bPart : particle name");
  brCmd->SetGuidance("  bProb : probability of Russian roulette");
  brCmd->SetGuidance("  bEnergy : max energy of a secondary for this biasing method");

  G4UIparameter* bPart = new G4UIparameter("bPart",'s',false);
  brCmd->SetParameter(bPart);

  G4UIparameter* bProb = new G4UIparameter("bProb",'d',false);
  brCmd->SetParameter(bProb);

  G4UIparameter* bEnergy = new G4UIparameter("bEnergy",'d',false);
  brCmd->SetParameter(bEnergy);

  G4UIparameter* bUnit = new G4UIparameter("bUnit",'s',true);
  brCmd->SetParameter(bUnit);
  brCmd->SetGuidance("unit of energy");

  brCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingMessenger::~StackingMessenger()
{
  delete killCmd;
  delete kCmd;
  delete brCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{     
  if(command == killCmd) {
    stackAction->SetKillStatus(killCmd->GetNewBoolValue(newValue));               
  } else if(command == kCmd) {
    stackAction->SetKill(newValue);

  } else if(command == brCmd) {
    G4double fb(1.0),en(1.e+30);
    G4String s1(""),unt("MeV");
    std::istringstream is(newValue);
    is >> s1 >> fb >> en >> unt;
    en *= G4UIcommand::ValueOf(unt);    
    stackAction->ActivateSecondaryBiasing(s1,fb,en);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
