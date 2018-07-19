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
/// \file PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class
//
// $Id: PrimaryGeneratorMessenger.cc 96158 2016-03-21 11:08:26Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                                   PrimaryGeneratorAction* Gun)
:G4UImessenger(),
  fAction(Gun),fGunDir(nullptr),fDefaultCmd(nullptr),
  fRndmCmd(nullptr),fTimeCmd(nullptr)
{ 
  fGunDir = new G4UIdirectory("/testhadr/gun/");
  fGunDir->SetGuidance("gun control");

  fDefaultCmd = new G4UIcmdWithoutParameter("/testhadr/gun/setDefault",this);
  fDefaultCmd->SetGuidance("set/reset kinematic defined in PrimaryGenerator");
  fDefaultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fRndmCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/gun/beamSize",this);
  fRndmCmd->SetGuidance("random lateral extension on the beam");
  fRndmCmd->SetParameterName("rBeam",false);
  fRndmCmd->SetRange("rBeam>=0.");
  fRndmCmd->SetUnitCategory("Length");
  fRndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fTimeCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/gun/beamTime",this);
  fTimeCmd->SetGuidance("set time of exposure");
  fTimeCmd->SetParameterName("tBeam",false);
  fTimeCmd->SetRange("tBeam>=0.");
  fTimeCmd->SetUnitCategory("Time");
  fTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fDefaultCmd;
  delete fRndmCmd;
  delete fTimeCmd;
  delete fGunDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if (command == fDefaultCmd)
   {fAction->SetDefaultKinematic();}

  if (command == fRndmCmd)
   {fAction->SetRndmBeam(fRndmCmd->GetNewDoubleValue(newValue));}

  if (command == fTimeCmd)
   {fAction->SetTimeExposure(fTimeCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
