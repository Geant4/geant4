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
/// \file electromagnetic/TestEm8/src/PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class
//
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
  fPhysDir = new G4UIdirectory("/testem/phys/");
  fPhysDir->SetGuidance("physics list commands");

  fECmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/setMaxE",this);  
  fECmd->SetGuidance("Set max energy deposit");
  fECmd->SetParameterName("Emax",false);
  fECmd->SetUnitCategory("Energy");
  fECmd->SetRange("Emax>0.0");
  fECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fEBCmd = new G4UIcmdWithAnInteger("/testem/phys/setNbinsE",this);  
  fEBCmd->SetGuidance("Set number of bins in energy.");
  fEBCmd->SetParameterName("Ebins",false);
  fEBCmd->SetRange("Ebins>0");
  fEBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fCBCmd = new G4UIcmdWithAnInteger("/testem/phys/setNbinsCl",this);  
  fCBCmd->SetGuidance("Set number of bins of clusters.");
  fCBCmd->SetParameterName("Cbins",false);
  fCBCmd->SetRange("Cbins>0");
  fCBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fCMCmd = new G4UIcmdWithAnInteger("/testem/phys/setMaxCl",this);  
  fCMCmd->SetGuidance("Set max number of clusters.");
  fCMCmd->SetParameterName("Cmax",false);
  fCMCmd->SetRange("Cmax>0");
  fCMCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fListCmd = new G4UIcmdWithAString("/testem/phys/addPhysics",this);  
  fListCmd->SetGuidance("Add modula physics list.");
  fListCmd->SetParameterName("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);

  fADCCmd = 
    new G4UIcmdWithADoubleAndUnit("/testem/phys/setEnergyPerChannel",this);
  fADCCmd->SetGuidance("Set energy per ADC channel");
  fADCCmd->SetParameterName("enadc",false,false);
  fADCCmd->SetUnitCategory("Energy");
  fADCCmd->SetDefaultUnit("keV");
  fADCCmd->SetRange("enadc>0.");
  fADCCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fNorCmd = new G4UIcmdWithADouble("/testem/phys/setNormFactor",this);
  fNorCmd->SetGuidance("Set factor for histogram normalisation");
  fNorCmd->SetParameterName("nfac",false,false);
  fNorCmd->SetRange("nfac>0.");
  fNorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSmCmd = new G4UIcmdWithADoubleAndUnit("/testem/phys/setEnergySmear",this);
  fSmCmd->SetGuidance("Set intrinsic width of detector response");
  fSmCmd->SetParameterName("sm",false,false);
  fSmCmd->SetDefaultUnit("keV");
  fSmCmd->SetRange("sm>0.");
  fSmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fECmd;
  delete fEBCmd;
  delete fCBCmd;
  delete fCMCmd;
  delete fListCmd;
  delete fADCCmd;
  delete fNorCmd;
  delete fSmCmd;
  delete fPhysDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  TestParameters* man = TestParameters::GetPointer();

  if( command == fECmd )
   { man->SetMaxEnergy(fECmd->GetNewDoubleValue(newValue)); }
  if( command == fEBCmd )
   { man->SetNumberBins(fEBCmd->GetNewIntValue(newValue)); }
  if( command == fCBCmd )
   { man->SetNumberBinsCluster(fCBCmd->GetNewIntValue(newValue)); }
  if( command == fCMCmd )
   { man->SetMaxCluster(fCMCmd->GetNewIntValue(newValue)); }
  if( command == fListCmd )
   { fPhysicsList->AddPhysicsList(newValue); }
  if( command == fADCCmd )
    { man->SetEnergyPerChannel(fADCCmd->GetNewDoubleValue(newValue)); }
  if( command == fNorCmd )
    { man->SetNormFactor(fNorCmd->GetNewDoubleValue(newValue)); }
  if( command == fSmCmd )
    { man->SetEnergySmear(fSmCmd->GetNewDoubleValue(newValue)); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
