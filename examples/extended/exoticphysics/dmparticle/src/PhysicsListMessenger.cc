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
/// \file exoticphysics/dmparticle/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
// $Id: PhysicsListMessenger.cc 95713 2016-02-22 08:08:38Z gcosmo $
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
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "TestParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:G4UImessenger(),fPhysicsList(pPhys)
{   
  fPhysDir = new G4UIdirectory("/testex/phys/");
  fPhysDir->SetGuidance("physics list commands");

  fECmd = new G4UIcmdWithADoubleAndUnit("/testex/phys/setMaxE",this);  
  fECmd->SetGuidance("Set max energy deposit");
  fECmd->SetParameterName("Emax",false);
  fECmd->SetUnitCategory("Energy");
  fECmd->SetRange("Emax>0.0");
  fECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fEBCmd = new G4UIcmdWithAnInteger("/testex/phys/setNbinsE",this);  
  fEBCmd->SetGuidance("Set number of bins in energy.");
  fEBCmd->SetParameterName("Ebins",false);
  fEBCmd->SetRange("Ebins>0");
  fEBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fListCmd = new G4UIcmdWithAString("/testex/phys/addPhysics",this);  
  fListCmd->SetGuidance("Add modula physics list.");
  fListCmd->SetParameterName("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);

  fPHCmd = new G4UIcmdWithADoubleAndUnit("/testex/phys/setLDMPhotonMass",this);
  fPHCmd->SetGuidance("Set LDMPhotonMass");
  fPHCmd->SetParameterName("LDMPmass",false);
  fPHCmd->SetUnitCategory("Energy");
  fPHCmd->SetRange("LDMPmass>0.0");
  fPHCmd->AvailableForStates(G4State_PreInit);

  fHCmd = new G4UIcmdWithADoubleAndUnit("/testex/phys/setLDMHiMass",this);
  fHCmd->SetGuidance("Set LDMHiMass");
  fHCmd->SetParameterName("LDMHmass",false);
  fHCmd->SetUnitCategory("Energy");
  fHCmd->SetRange("LDMHmass>0.0");
  fHCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fECmd;
  delete fPHCmd;
  delete fHCmd;
  delete fEBCmd;
  delete fListCmd;
  delete fPhysDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                       G4String newValue)
{       
  TestParameters* man = TestParameters::GetPointer();

  if( command == fECmd )
  { 
    man->SetMaxEnergy(fECmd->GetNewDoubleValue(newValue)); 
  }
  if( command == fEBCmd )
  { 
    man->SetNumberBins(fEBCmd->GetNewIntValue(newValue)); 
  }
  if( command == fListCmd )
  { 
    fPhysicsList->AddPhysicsList(newValue); 
  }
  if( command == fPHCmd )
  { 
    fPhysicsList->SetLDMPhotonMass(fPHCmd->GetNewDoubleValue(newValue)); 
  }
  if( command == fHCmd )
  { 
    fPhysicsList->SetLDMHiMass(fHCmd->GetNewDoubleValue(newValue)); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

