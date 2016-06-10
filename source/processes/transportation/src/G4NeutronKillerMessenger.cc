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
// $Id: G4NeutronKillerMessenger.cc 68048 2013-03-13 14:34:07Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4NeutronKillerMessenger
//
// Description: Messenger class
//
// Author:      V.Ivanchenko 26/09/00 for HARP software
//
//----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4NeutronKillerMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIdirectory.hh"
#include "G4NeutronKiller.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NeutronKillerMessenger::G4NeutronKillerMessenger(G4NeutronKiller* p)
  :killer(p)
{

  dir = new G4UIdirectory("/physics_engine/neutron/");
  dir->SetGuidance("control on neutrons");

  eCmd = new G4UIcmdWithADoubleAndUnit("/physics_engine/neutron/energyLimit",this);
  eCmd->SetGuidance("Set tracking cut - min energy of a particle.");
  eCmd->SetParameterName("energyLimit",false);
  eCmd->SetUnitCategory("Energy");
  eCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  tCmd = new G4UIcmdWithADoubleAndUnit("/physics_engine/neutron/timeLimit",this);
  tCmd->SetGuidance("Set time limit.");
  tCmd->SetParameterName("timeLimit",false);
  tCmd->SetUnitCategory("Time");
  tCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NeutronKillerMessenger::~G4NeutronKillerMessenger()
{
  delete eCmd;
  delete tCmd;
  delete dir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NeutronKillerMessenger::SetNewValue(G4UIcommand* command, G4String val)
{ 
  if (command == eCmd)
    killer->SetKinEnergyLimit(eCmd->GetNewDoubleValue(val));

  if (command == tCmd)
    killer->SetTimeLimit(tCmd->GetNewDoubleValue(val));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
