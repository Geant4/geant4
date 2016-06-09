//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: PhysicsList.cc,v 1.7 2005/06/07 13:55:07 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Author:      V.Ivanchenko 03.05.2004
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "ParticlesBuilder.hh"
#include "G4EmQEDBuilder.hh"
#include "G4EmMuonBuilder.hh"
#include "G4EmHadronBuilder.hh"
#include "G4EmQEDBuilder52.hh"
#include "G4LowEnergyQEDBuilder.hh"
#include "G4PenelopeQEDBuilder.hh"
#include "G4EmMuonBuilder52.hh"
#include "G4EmHadronBuilder52.hh"
#include "G4StepLimiterBuilder.hh"
#include "DecaysBuilder.hh"
#include "EmHadronElasticBuilder.hh"
#include "EmBinaryCascadeBuilder.hh"
#include "EmIonBinaryCascadeBuilder.hh"
#include "EmGammaNucleusBuilder.hh"

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
  verbose = 0;
  G4LossTableManager::Instance()->SetVerbose(0);
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  pMessenger = new PhysicsListMessenger(this);

  // Add Physics builders
  RegisterPhysics(new ParticlesBuilder());
  steplimiter = new G4StepLimiterBuilder();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(verbose > 0) 
    G4cout << "Construte Particles" << G4endl;
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  if(verbose > 0) 
    G4cout << "Construte Processes" << G4endl;
  if(!emBuilderIsRegisted) AddPhysicsList("standard");
  G4VModularPhysicsList::ConstructProcess();

  // Define energy interval for loss processes
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(0.1*keV);
  emOptions.SetMaxEnergy(100.*GeV);
  emOptions.SetDEDXBinning(90);
  emOptions.SetLambdaBinning(90);
  emOptions.SetBuildPreciseRange(false);
  emOptions.SetApplyCuts(true);
  //emOptions.SetVerbose(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if(verbose > 0) {
    G4cout << "Add Physics <" << name 
           << "> emBuilderIsRegisted= " << emBuilderIsRegisted
           << G4endl;
  }
  if ((name == "standard") && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmQEDBuilder());
    RegisterPhysics(steplimiter);
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmHadronBuilder());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;    

  } else if (name == "g4v52" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmQEDBuilder52());
    RegisterPhysics(steplimiter);
    RegisterPhysics(new G4EmMuonBuilder52());
    RegisterPhysics(new G4EmHadronBuilder52());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "lowenergy" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4LowEnergyQEDBuilder());
    RegisterPhysics(steplimiter);
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmHadronBuilder());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "penelope" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4PenelopeQEDBuilder());
    RegisterPhysics(steplimiter);
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmHadronBuilder());
    emBuilderIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "decay" && !decayIsRegisted) {
    RegisterPhysics(new DecaysBuilder());
    decayIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "elastic" && !helIsRegisted) {
    RegisterPhysics(new EmHadronElasticBuilder());
    helIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    
  } else if (name == "binary" && !bicIsRegisted) {
    RegisterPhysics(new EmBinaryCascadeBuilder());
    bicIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    
  } else if (name == "binary_ion" && !ionIsRegisted) {
    RegisterPhysics(new EmIonBinaryCascadeBuilder());
    ionIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;

  } else if (name == "gamma_nuc" && !gnucIsRegisted) {
    RegisterPhysics(new EmGammaNucleusBuilder());
    gnucIsRegisted = true;
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
