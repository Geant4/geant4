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
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "OpNovicePhysicsListMessenger.hh"

#include "OpNovicePhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNovicePhysicsListMessenger::
  OpNovicePhysicsListMessenger(OpNovicePhysicsList* pPhys):pPhysicsList(pPhys)
{
  OpNoviceDir = new G4UIdirectory("/OpNovice/");
  OpNoviceDir->SetGuidance("UI commands of this example");
  
  physDir = new G4UIdirectory("/OpNovice/phys/");
  physDir->SetGuidance("PhysicsList control");
 
  verboseCmd = new G4UIcmdWithAnInteger("/OpNovice/phys/verbose",this);  
  verboseCmd->SetGuidance("set verbose for physics processes");
  verboseCmd->SetParameterName("verbose",true);
  verboseCmd->SetDefaultValue(1);
  verboseCmd->SetRange("verbose>=0");
  verboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  cerenkovCmd =
           new G4UIcmdWithAnInteger("/OpNovice/phys/cerenkovMaxPhotons",this);  
  cerenkovCmd->SetGuidance("set max nb of photons per step");
  cerenkovCmd->SetParameterName("MaxNumber",false);
  cerenkovCmd->SetRange("MaxNumber>=0");
  cerenkovCmd->AvailableForStates(G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNovicePhysicsListMessenger::~OpNovicePhysicsListMessenger()
{
  delete verboseCmd;
  delete cerenkovCmd;
  delete physDir;
  delete OpNoviceDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNovicePhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  if( command == verboseCmd )
   { pPhysicsList->SetVerbose(verboseCmd->GetNewIntValue(newValue));}
   
  if( command == cerenkovCmd )
   {pPhysicsList->SetNbOfPhotonsCerenkov(cerenkovCmd->GetNewIntValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
