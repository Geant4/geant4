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


#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:fPL(pPhys)
{
  fPhysDir = new G4UIdirectory("/my_phys/");
  fPhysDir->SetGuidance("Commands to set physics list");

  fPhysListCmd = new G4UIcmdWithAString("/my_phys/setList", this);
  fPhysListCmd->SetGuidance("Set a *reference* physics list.");
  fPhysListCmd->SetParameterName("physList", false);
  fPhysListCmd->AvailableForStates(G4State_PreInit);

  fVerbCmd = new G4UIcmdWithAnInteger("/my_phys/verbose", this);
  fVerbCmd->SetGuidance("Set verbose level for physics list");
  fVerbCmd->SetParameterName("verb", false);
  fVerbCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fPhysListCmd;
  delete fVerbCmd;
  delete fPhysDir;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if ( command == fPhysListCmd )
    fPL->SetPhysicsList(newValue);

  else if ( command == fVerbCmd )
    fPL->SetVerbose(fVerbCmd->GetNewIntValue(newValue));
}

