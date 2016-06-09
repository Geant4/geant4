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
// $Id: PhysicsListMessenger.cc,v 1.1 2007/10/15 16:20:23 maire Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:pPhysicsList(pPhys)
{
  physDir = new G4UIdirectory("/testem/phys/");
  physDir->SetGuidance("physics list commands");
  
  singleScat = new G4UIcmdWithABool("/testem/phys/singleScattering",this);
  singleScat->SetGuidance("apply single Coulomb scattering process");
  singleScat->SetParameterName("flag",true);
  singleScat->SetDefaultValue(true);
  singleScat->AvailableForStates(G4State_PreInit);
        
  brem = new G4UIcmdWithABool("/testem/phys/registerBrem",this);
  brem->SetGuidance("register Brems in PhysicsList");
  brem->SetParameterName("flag",true);
  brem->SetDefaultValue(true);
  brem->AvailableForStates(G4State_PreInit);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete singleScat;
  delete brem;
  delete physDir;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{  
  if (command == singleScat)
   {pPhysicsList->SingleCoulombScattering(singleScat->GetNewBoolValue(newValue));}
            
  if (command == brem)
   {pPhysicsList->RegisterBrem(brem->GetNewBoolValue(newValue));}      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
