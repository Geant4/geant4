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
#include "G4EmStandardPhysics_option4.hh"
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
#include "G4EmParameters.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
  fHelIsRegisted = false;
  fBicIsRegisted = false;
  fIonIsRegisted = false;
  fGnucIsRegisted = false;
  fStopIsRegisted = false;
  fVerbose = 1;

  SetDefaultCutValue(1*mm);

  fMessenger = new PhysicsListMessenger(this);

  // Add Physics builders
  RegisterPhysics(new G4EmStandardPhysics());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new StepLimiterBuilder());
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

  G4VModularPhysicsList::ConstructProcess();

  // Define energy interval for loss processes
  // from 10 eV to 10 GeV
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetMinEnergy(0.01*keV);
  param->SetMaxEnergy(10.*GeV);
  param->SetNumberOfBinsPerDecade(10);
  param->SetVerbose(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if(fVerbose > 0) {
    G4cout << "### PhysicsList Add Physics <" << name 
           << "> " << G4endl;
  }
  if (name == "emstandard") {
    ReplacePhysics(new G4EmStandardPhysics());

  } else if (name == "emstandard_opt1") {
    ReplacePhysics(new G4EmStandardPhysics_option1());

  } else if (name == "emstandard_opt2") {
    ReplacePhysics(new G4EmStandardPhysics_option2());

  } else if (name == "emstandard_opt3") {
    ReplacePhysics(new G4EmStandardPhysics_option3());

  } else if (name == "emstandard_opt4") {
    ReplacePhysics(new G4EmStandardPhysics_option4());

  } else if (name == "emlivermore") {
    ReplacePhysics(new G4EmLivermorePhysics());

  } else if (name == "empenelope") {
    ReplacePhysics(new G4EmPenelopePhysics());

  } else if (name == "emlowenergy") {
    ReplacePhysics(new G4EmLowEPPhysics());

  } else if (name == "elastic" && !fHelIsRegisted) {
    RegisterPhysics(new G4HadronElasticPhysics());
    fHelIsRegisted = true;
    
  } else if (name == "binary" && !fBicIsRegisted) {
    RegisterPhysics(new G4HadronInelasticQBBC());
    fBicIsRegisted = true;
    
  } else if (name == "binary_ion" && !fIonIsRegisted) {
    RegisterPhysics(new G4IonBinaryCascadePhysics());
    fIonIsRegisted = true;

  } else if (name == "gamma_nuc" && !fGnucIsRegisted) {
    RegisterPhysics(new G4EmExtraPhysics());
    fGnucIsRegisted = true;

  } else if (name == "stopping" && !fStopIsRegisted) {
    RegisterPhysics(new G4StoppingPhysics());
    fStopIsRegisted = true;
    
  } else {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" 
           << " fail - module is already regitered or is unknown " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
