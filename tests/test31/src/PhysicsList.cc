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
// $Id: PhysicsList.cc,v 1.11 2004-11-03 12:36:46 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4EmQEDBuilder70.hh"
#include "G4EmLowEnergyQEDBuilder.hh"
#include "G4EmPenelopeQEDBuilder.hh"
#include "G4EmMuonBuilder.hh"
#include "G4EmHadronBuilder.hh"
#include "G4EmLowEnergyHadronBuilder.hh"
#include "G4EmLowEnergyHadronBuilderMA.hh"
#include "G4EmHighEnergyBuilder.hh"
#include "G4EmQEDBuilder52.hh"
#include "G4EmMuonBuilder52.hh"
#include "G4EmHadronBuilder52.hh"
#include "DecaysBuilder.hh"

#include "PhysListBinaryCascade.hh"
#include "PhysListHadronElastic.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList(),
  verbose(0),
  emBuilderIsRegisted(false)
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  pMessenger = new PhysicsListMessenger(this);

  // Add Physics builders
  RegisterPhysics(new ParticlesBuilder());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(!emBuilderIsRegisted) AddPhysicsList("standard");
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();

  // Define energy interval for loss processes
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(0.1*keV);
  //  emOptions.SetMaxEnergy(100.*GeV);
  emOptions.SetMaxEnergy(1.*TeV);
  emOptions.SetDEDXBinning(100);
  emOptions.SetLambdaBinning(100);
  emOptions.SetBuildPreciseRange(false);
  emOptions.SetVerbose(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  G4bool yes = false;
  if (verbose>-1) {
    G4cout << "PhysicsList::Add New PhysicsList <" << name << ">  " 
           <<  emBuilderIsRegisted<< G4endl;
  }

  if ("standard" == name && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmQEDBuilder());
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmHadronBuilder());
    RegisterPhysics(new DecaysBuilder());
    emBuilderIsRegisted = true;
    yes = true;

  } else if ("standard70" == name  && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmQEDBuilder70());
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmHadronBuilder());
    RegisterPhysics(new DecaysBuilder());
    emBuilderIsRegisted = true;
    yes = true;

  } else if ("g4v52" == name && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmQEDBuilder52());
    RegisterPhysics(new G4EmMuonBuilder52());
    RegisterPhysics(new G4EmHadronBuilder52());
    RegisterPhysics(new DecaysBuilder());
    emBuilderIsRegisted = true;
    yes = true;

  } else if ("lowenergy" == name && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmLowEnergyQEDBuilder());
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmLowEnergyHadronBuilder());
    RegisterPhysics(new DecaysBuilder());
    emBuilderIsRegisted = true;
    yes = true;

  } else if ("lowenergyMA" == name && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmQEDBuilder());
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmLowEnergyHadronBuilderMA());
    RegisterPhysics(new DecaysBuilder());
    emBuilderIsRegisted = true;
    yes = true;

  } else if ("penelope" == name && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmPenelopeQEDBuilder());
    RegisterPhysics(new G4EmMuonBuilder());
    RegisterPhysics(new G4EmLowEnergyHadronBuilder());
    RegisterPhysics(new DecaysBuilder());
    emBuilderIsRegisted = true;
    yes = true;

  } else if ("high_energy" == name ) {
    RegisterPhysics(new G4EmHighEnergyBuilder());
    yes = true;

  } else if ("binary" == name ) {
    RegisterPhysics(new PhysListBinaryCascade());
    yes = true;

  } else if ("elastic" == name ) {
    RegisterPhysics(new PhysListHadronElastic());
    yes = true;

  } else {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">"
           << " fail - name is unknown " << G4endl;
  }
  if (verbose>-1 && yes) {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" << G4endl;
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
