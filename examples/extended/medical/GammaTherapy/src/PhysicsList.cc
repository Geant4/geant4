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
// $Id: PhysicsList.cc,v 1.17 2009-11-21 16:47:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4StepLimiterBuilder.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4QStoppingPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
: G4VModularPhysicsList()
{
  emBuilderIsRegisted = false;
  decayIsRegisted = false;
  stepLimiterIsRegisted = false;
  helIsRegisted = false;
  bicIsRegisted = false;
  ionIsRegisted = false;
  gnucIsRegisted = false;
  stopIsRegisted = false;
  verbose = 1;
  G4LossTableManager::Instance()->SetVerbose(1);
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  pMessenger = new PhysicsListMessenger(this);

  // Add Physics builders
  RegisterPhysics(new G4DecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(verbose > 0) 
    G4cout << "### PhysicsList Construte Particles" << G4endl;
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  if(verbose > 0) 
    G4cout << "### PhysicsList Construte Processes" << G4endl;
  if(!emBuilderIsRegisted) AddPhysicsList("emstandard");
  RegisterPhysics(new G4StepLimiterBuilder());
  G4VModularPhysicsList::ConstructProcess();

  // Define energy interval for loss processes
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(0.1*keV);
  emOptions.SetMaxEnergy(100.*GeV);
  emOptions.SetDEDXBinning(90);
  emOptions.SetLambdaBinning(90);
  //  emOptions.SetBuildCSDARange(false);
  emOptions.SetApplyCuts(true);
  //emOptions.SetVerbose(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if(verbose > 0) {
    G4cout << "### PhysicsList Add Physics <" << name 
           << "> emBuilderIsRegisted= " << emBuilderIsRegisted
           << G4endl;
  }
  if ((name == "emstandard") && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;    

  } else if (name == "emstandard_opt1" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option1());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "emstandard_opt2" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option2());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "emstandard_opt3" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option3());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "emlivermore" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmLivermorePhysics());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "empenelope" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmPenelopePhysics());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "elastic" && !helIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4HadronElasticPhysics());
    helIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    
  } else if (name == "binary" && !bicIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4HadronInelasticQBBC());
    bicIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    
  } else if (name == "binary_ion" && !ionIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4IonBinaryCascadePhysics());
    ionIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "gamma_nuc" && !gnucIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4EmExtraPhysics());
    gnucIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "stopping" && !stopIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4QStoppingPhysics());
    gnucIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    
  } else if(!emBuilderIsRegisted) {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" 
           << " fail - EM physics should be registered first " << G4endl;
  } else {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" 
           << " fail - module is already regitered or is unknown " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  if (verbose>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
