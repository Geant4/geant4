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
/// \file medical/GammaTherapy/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 82277 2014-06-13 14:40:54Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Author:      V.Ivanchenko 03.05.2004
//
// Modified:
// 16.11.06 Use components from physics_lists subdirectory (V.Ivanchenko)
// 16.05.07 Use renamed EM components from physics_lists (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "StepLimiterBuilder.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
  fEmBuilderIsRegisted = false;
  fHelIsRegisted = false;
  fBicIsRegisted = false;
  fIonIsRegisted = false;
  fGnucIsRegisted = false;
  fStopIsRegisted = false;
  fVerbose = 1;
  G4LossTableManager::Instance()->SetVerbose(fVerbose);
  SetDefaultCutValue(1*mm);

  fMessenger = new PhysicsListMessenger(this);

  // Add Physics builders
  RegisterPhysics(new G4DecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(fVerbose > 0) {
    G4cout << "### PhysicsList Construte Particles" << G4endl;
  }
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  if(fVerbose > 0) {
    G4cout << "### PhysicsList Construte Processes" << G4endl;
  }
  if(!fEmBuilderIsRegisted) { AddPhysicsList("emstandard"); }
  RegisterPhysics(new StepLimiterBuilder());
  G4VModularPhysicsList::ConstructProcess();

  // Define energy interval for loss processes
  // from 10 eV to 10 GeV
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(0.01*keV);
  emOptions.SetMaxEnergy(10.*GeV);
  emOptions.SetDEDXBinning(90);
  emOptions.SetLambdaBinning(90);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if(fVerbose > 0) {
    G4cout << "### PhysicsList Add Physics <" << name 
           << "> " << G4endl;
  }
  if ((name == "emstandard") && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics());
    fEmBuilderIsRegisted = true;

  } else if (name == "emstandard_opt1" && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option1());
    fEmBuilderIsRegisted = true;

  } else if (name == "emstandard_opt2" && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option2());
    fEmBuilderIsRegisted = true;

  } else if (name == "emstandard_opt3" && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option3());
    fEmBuilderIsRegisted = true;

  } else if (name == "emlivermore" && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmLivermorePhysics());
    fEmBuilderIsRegisted = true;

  } else if (name == "empenelope" && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmPenelopePhysics());
    fEmBuilderIsRegisted = true;

  } else if (name == "emlowenergy" && !fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmLowEPPhysics());
    fEmBuilderIsRegisted = true;

  } else if (name == "elastic" && !fHelIsRegisted && fEmBuilderIsRegisted) {
    RegisterPhysics(new G4HadronElasticPhysics());
    fHelIsRegisted = true;
    
  } else if (name == "binary" && !fBicIsRegisted && fEmBuilderIsRegisted) {
    RegisterPhysics(new G4HadronInelasticQBBC());
    fBicIsRegisted = true;
    
  } else if (name == "binary_ion" && !fIonIsRegisted && fEmBuilderIsRegisted) {
    RegisterPhysics(new G4IonBinaryCascadePhysics());
    fIonIsRegisted = true;

  } else if (name == "gamma_nuc" && !fGnucIsRegisted && fEmBuilderIsRegisted) {
    RegisterPhysics(new G4EmExtraPhysics());
    fGnucIsRegisted = true;

  } else if (name == "stopping" && !fStopIsRegisted && fEmBuilderIsRegisted) {
    RegisterPhysics(new G4StoppingPhysics());
    fStopIsRegisted = true;
    
  } else if(!fEmBuilderIsRegisted) {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" 
           << " fail - EM physics should be registered first " << G4endl;
  } else {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" 
           << " fail - module is already regitered or is unknown " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
