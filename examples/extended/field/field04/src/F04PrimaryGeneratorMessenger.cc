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
/// \file field/field04/src/F04PrimaryGeneratorMessenger.cc
/// \brief Implementation of the F04PrimaryGeneratorMessenger class
//

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4SystemOfUnits.hh"

#include "F04PrimaryGeneratorAction.hh"
#include "F04PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PrimaryGeneratorMessenger::
                  F04PrimaryGeneratorMessenger(F04PrimaryGeneratorAction* gun)
  : fAction(gun)
{
  fRndmCmd = new G4UIcmdWithAString("/gun/random",this);
  fRndmCmd->SetGuidance("Shoot randomly the incident particle.");
  fRndmCmd->SetGuidance("  Choice : on, off(default)");
  fRndmCmd->SetParameterName("choice",true);
  fRndmCmd->SetDefaultValue("off");
  fRndmCmd->SetCandidates("on off");
  fRndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fSetXvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/xvertex",this);
  fSetXvertexCmd->SetGuidance(" Set x coord. of the primary vertex.");
  fSetXvertexCmd->SetParameterName("xv",true);
  fSetXvertexCmd->SetDefaultValue(0.0*mm);
  fSetXvertexCmd->SetDefaultUnit("mm");
 
  fSetYvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/yvertex",this);
  fSetYvertexCmd->SetGuidance(" Set y coord. of the primary vertex.");
  fSetYvertexCmd->SetParameterName("yv",true);
  fSetYvertexCmd->SetDefaultValue(0.0*mm);
  fSetYvertexCmd->SetDefaultUnit("mm");
 
  fSetZvertexCmd = new G4UIcmdWithADoubleAndUnit("/gun/zvertex",this);
  fSetZvertexCmd->SetGuidance(" Set z coord. of the primary vertex.");
  fSetZvertexCmd->SetParameterName("zv",true);
  fSetZvertexCmd->SetDefaultValue(0.0*mm);
  fSetZvertexCmd->SetDefaultUnit("mm");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PrimaryGeneratorMessenger::~F04PrimaryGeneratorMessenger()
{
  delete fRndmCmd;
  delete fSetXvertexCmd;
  delete fSetYvertexCmd;
  delete fSetZvertexCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PrimaryGeneratorMessenger::
                          SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command == fRndmCmd )
   { fAction->SetRndmFlag(newValue);}
  if( command == fSetXvertexCmd)
   { fAction->SetXvertex(fSetXvertexCmd->GetNewDoubleValue(newValue));}
  if( command == fSetYvertexCmd)
   { fAction->SetYvertex(fSetYvertexCmd->GetNewDoubleValue(newValue));}
  if( command == fSetZvertexCmd)
   { fAction->SetZvertex(fSetZvertexCmd->GetNewDoubleValue(newValue));}
}
