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
/// \file electromagnetic/TestEm5/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsListMessenger.hh"
#include "G4SystemOfUnits.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
  :G4UImessenger(),fPhysicsList(pPhys),fMaxChargedStep(1*CLHEP::mm)
{
  fPhysDir = new G4UIdirectory("/testem/phys/");
  fPhysDir->SetGuidance("physics list commands");
<<<<<<< HEAD
   
  fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/microbeam/phys/setGCut",this);  
  fGammaCutCmd->SetGuidance("Set gamma cut.");
  fGammaCutCmd->SetParameterName("Gcut",false);
  fGammaCutCmd->SetUnitCategory("Length");
  fGammaCutCmd->SetRange("Gcut>0.0");
  fGammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fElectCutCmd = new G4UIcmdWithADoubleAndUnit("/microbeam/phys/setECut",this);  
  fElectCutCmd->SetGuidance("Set electron cut.");
  fElectCutCmd->SetParameterName("Ecut",false);
  fElectCutCmd->SetUnitCategory("Length");
  fElectCutCmd->SetRange("Ecut>0.0");
  fElectCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fProtoCutCmd = new G4UIcmdWithADoubleAndUnit("/microbeam/phys/setPCut",this);  
  fProtoCutCmd->SetGuidance("Set positron cut.");
  fProtoCutCmd->SetParameterName("Pcut",false);
  fProtoCutCmd->SetUnitCategory("Length");
  fProtoCutCmd->SetRange("Pcut>0.0");
  fProtoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fAllCutCmd = new G4UIcmdWithADoubleAndUnit("/microbeam/phys/setCuts",this);  
  fAllCutCmd->SetGuidance("Set cut for all.");
  fAllCutCmd->SetParameterName("cut",false);
  fAllCutCmd->SetUnitCategory("Length");
  fAllCutCmd->SetRange("cut>0.0");
  fAllCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fListCmd = new G4UIcmdWithAString("/microbeam/phys/addPhysics",this);  
=======

  fListCmd = new G4UIcmdWithAString("/testem/phys/addPhysics",this);  
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  fListCmd->SetGuidance("Add modula physics list.");
  fListCmd->SetParameterName("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);
  fListCmd->SetToBeBroadcasted(false);      

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/testem/stepMax",this);
  fStepMaxCmd->SetGuidance("Set max allowed step length");
  fStepMaxCmd->SetParameterName("mxStep",false);
  fStepMaxCmd->SetRange("mxStep>0.");
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fGammaCutCmd;
  delete fElectCutCmd;
  delete fProtoCutCmd;
  delete fAllCutCmd;
  delete fListCmd;
  delete fPhysDir;
  delete fStepMaxCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

<<<<<<< HEAD
void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  if( command == fGammaCutCmd )
   { fPhysicsList->SetCutForGamma(fGammaCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == fElectCutCmd )
   { fPhysicsList->SetCutForElectron(fElectCutCmd->GetNewDoubleValue(newValue));}
     
  if( command == fProtoCutCmd )
   { fPhysicsList->SetCutForPositron(fProtoCutCmd->GetNewDoubleValue(newValue));}

  if( command == fAllCutCmd )
   {
     G4double cut = fAllCutCmd->GetNewDoubleValue(newValue);
     fPhysicsList->SetCutForGamma(cut);
     fPhysicsList->SetCutForElectron(cut);
     fPhysicsList->SetCutForPositron(cut);
   } 

=======
void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
  if( command == fListCmd )
    { fPhysicsList->AddPhysicsList(newValue); }
  if (command == fStepMaxCmd)
    { fMaxChargedStep = fStepMaxCmd->GetNewDoubleValue(newValue); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
