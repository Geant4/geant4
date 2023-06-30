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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"

#include "PhysListEmStandard.hh"
#include "PhysicsListMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  fHelIsRegisted = false;
  fBicIsRegisted = false;
  fBiciIsRegisted = false;

  fMessenger = new PhysicsListMessenger(this);
  SetVerboseLevel(1);

  // EM physics
  fEmName = G4String("emlivermore");
  fEmPhysics = new G4EmLivermorePhysics();
  //  fEmPhysicsList = new PhysListEmStandard();
  //  fEmName = G4String("local");

  // Decay physics
  fDecayPhysics = new G4DecayPhysics(1);

  //  SetDefaultCutValue(1*mm);
  SetDefaultCutValue(0.01 * um);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fEmPhysics;
  delete fDecayPhysics;
  for (size_t i = 0; i < fHadronPhysics.size(); i++) {
    delete fHadronPhysics[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fDecayPhysics->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  // transportation
  //
  AddTransportation();

  // electromagnetic physics
  //
  fEmPhysics->ConstructProcess();

  // decay physics
  //
  fDecayPhysics->ConstructProcess();

  // hadronic physics
  for (size_t i = 0; i < fHadronPhysics.size(); i++) {
    fHadronPhysics[i]->ConstructProcess();
  }

  // step limitation (as a full process)
  //
  AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  auto fStepMaxProcess = new G4StepLimiter();

  auto particleIterator = GetParticleIterator();

  particleIterator->reset();

  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (fStepMaxProcess->IsApplicable(*particle)) {
      pmanager->AddDiscreteProcess(fStepMaxProcess);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel > -1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "local") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new PhysListEmStandard(name);
  }
  else if (name == "emstandard_opt0") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysics();
  }
  else if (name == "emstandard_opt1") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysics_option1();
  }
  else if (name == "emstandard_opt2") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysics_option2();
  }
  else if (name == "emstandard_opt3") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysics_option3();
  }
  else if (name == "emstandard_opt4") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysics_option4();
  }
  else if (name == "emstandardSS") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysicsSS();
  }
  else if (name == "emstandardWVI") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysicsWVI();
  }
  else if (name == "emstandardGS") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmStandardPhysicsGS(verboseLevel);
  }
  else if (name == "empenelope") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmPenelopePhysics();
  }
  else if (name == "emlowenergy") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmLowEPPhysics();
  }
  else if (name == "emlivermore") {
    fEmName = name;
    delete fEmPhysics;
    fEmPhysics = new G4EmLivermorePhysics();
  }
  else if (name == "elastic" && !fHelIsRegisted) {
    fHadronPhysics.push_back(new G4HadronElasticPhysics(verboseLevel));
    fHelIsRegisted = true;
  }
  else if (name == "DElastic" && !fHelIsRegisted) {
    fHadronPhysics.push_back(new G4HadronDElasticPhysics(verboseLevel));
    fHelIsRegisted = true;
  }
  else if (name == "HElastic" && !fHelIsRegisted) {
    fHadronPhysics.push_back(new G4HadronHElasticPhysics(verboseLevel));
    fHelIsRegisted = true;
  }
  else if (name == "binary" && !fBicIsRegisted) {
    fHadronPhysics.push_back(new G4HadronInelasticQBBC(verboseLevel));
    fBicIsRegisted = true;
  }
  else if (name == "binary_ion" && !fBiciIsRegisted) {
    fHadronPhysics.push_back(new G4IonPhysics(verboseLevel));
    fBiciIsRegisted = true;
  }
  else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
