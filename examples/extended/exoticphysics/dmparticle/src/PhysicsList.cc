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
/// \file exoticphysics/dmparticle/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.hh 92047 2015-08-14 07:23:37Z gcosmo $
//
/////////////////////////////////////////////////
//
// ClassName:   PhysicsList
//
// Authors:  01.06.17 V.Ivanchenko 
//
//
///////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4DecayPhysics.hh"

#include "G4LDMHi.hh"
#include "G4LDMHiBar.hh"
#include "G4LDMPhoton.hh"
#include "G4LDMBremsstrahlung.hh"
#include "G4LDMBremModel.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmConfigurator.hh"
#include "G4EmParameters.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
  :G4VModularPhysicsList(),
   fLDMPhotonMass(0.5*CLHEP::GeV),fLDMHiMass(0.1*CLHEP::GeV),
   fLDMPhoton(true),fLDMHi(false)
{ 
  fMessenger = new PhysicsListMessenger(this);

  // Decay Physics is always defined
  fDecayPhysicsList = new G4DecayPhysics();

  // EM physics
  fEmName = G4String("emstandard_opt0");
  fEmPhysicsList = new G4EmStandardPhysics(1);

  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fDecayPhysicsList;
  delete fEmPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  SetDefaultCutValue(1*mm);
  fDecayPhysicsList->ConstructParticle();
  if(fLDMPhoton) { G4LDMPhoton::LDMPhotonDefinition(fLDMPhotonMass); }
  if(fLDMHi) { 
    G4LDMHi::LDMHiDefinition(fLDMHiMass);
    G4LDMHiBar::LDMHiBarDefinition(fLDMHiMass);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess(); 
  fDecayPhysicsList->ConstructProcess();
  AddDarkMatter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  //G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  if (name == fEmName) { 
    return; 

  } else if (name == "emstandard_opt0") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();

  } else if (name == "emstandard_opt1") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "emstandard_opt2") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt3") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();

  } else if (name == "emstandard_opt4") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();

  } else if (name == "emstandardWVI") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsWVI();

  } else if (name == "emstandardSS") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsSS();

  } else if (name == "emstandardGS") {
    fEmName = name;
    delete fEmPhysicsList;
    
    fEmPhysicsList = new G4EmStandardPhysicsGS();

  } else if (name == "emlivermore") {
    fEmName = name;    
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();

  } else if (name == "empenelope") {
    fEmName = name;
    delete fEmPhysicsList;    
    fEmPhysicsList = new G4EmPenelopePhysics();

  } else if (name == "emlowenergy") {
    fEmName = name;
    delete fEmPhysicsList;    
    fEmPhysicsList = new G4EmLowEPPhysics();

  } else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100.*eV,1e5);
  if ( verboseLevel > 0 ) { DumpCutValuesTable(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddDarkMatter()
{
  G4cout << " PhysicsList::AddDarkMatter: " << fLDMPhoton << G4endl;
  if(fLDMPhoton) {
    //G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    G4LDMBremsstrahlung* ldmb = new G4LDMBremsstrahlung();
    G4ParticleDefinition* p = G4Proton::Proton();
    G4ProcessManager* man = p->GetProcessManager();
    man->AddProcess(ldmb, -1, -1, 5);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetLDMPhotonMass(G4double val) 
{ 
  G4cout << "### PhysicsList::SetLDMPhotonMass: new value " << val/GeV
         << " GeV" << G4endl;
  fLDMPhotonMass = val; 
  fLDMPhoton = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetLDMHiMass(G4double val)
{ 
  fLDMHiMass = val; 
  fLDMHi = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
