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
// $Id: PhysicsListMessenger.cc,v 1.3 2010-09-08 09:12:10 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsListMessenger
//
// Description: EM physics with a possibility to add PAI model
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:pPhysicsList(pPhys)
{   
  eCmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/setMaxE",this);  
  eCmd->SetGuidance("Set max energy deposit");
  eCmd->SetParameterName("Emax",false);
  eCmd->SetUnitCategory("Energy");
  eCmd->SetRange("Emax>0.0");
  eCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ebCmd = new G4UIcmdWithAnInteger("/testem/phys/setNbinsE",this);  
  ebCmd->SetGuidance("Set number of bins in energy.");
  ebCmd->SetParameterName("Ebins",false);
  ebCmd->SetRange("Ebins>0");
  ebCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cbCmd = new G4UIcmdWithAnInteger("/testem/phys/setNbinsCl",this);  
  cbCmd->SetGuidance("Set max number of clusters.");
  cbCmd->SetParameterName("Cbins",false);
  cbCmd->SetRange("Cbins>0");
  cbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  pListCmd = new G4UIcmdWithAString("/testem/phys/addPhysics",this);  
  pListCmd->SetGuidance("Add modula physics list.");
  pListCmd->SetParameterName("PList",false);
  pListCmd->AvailableForStates(G4State_PreInit);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete eCmd;
  delete ebCmd;
  delete cbCmd;
  delete pListCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  HistoManager* man = HistoManager::GetPointer();
  if( command == eCmd )
   { man->SetMaxEnergy(eCmd->GetNewDoubleValue(newValue));}
     
  if( command == ebCmd )
   { man->SetNumberBins(ebCmd->GetNewIntValue(newValue));}
     
  if( command == cbCmd )
   { man->SetNumberBinsCluster(cbCmd->GetNewIntValue(newValue));}

  if( command == pListCmd )
   { pPhysicsList->AddPhysicsList(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
