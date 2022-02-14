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
/// \file parallel/ThreadsafeScorers/src/TSPhysicsList.cc
/// \brief Implementation of the TSPhysicsList class
//
//
//
//
/// This is a very, very extensive physics list and step-limiters are applied
///     to many particles. The reasoning behind this is because we wan't to put
///     as much pressure on the atomics as possible and produce as much
///     round-off error as possible. See descriptions in README and
///     TSDetectorConstruction for more details.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TSPhysicsList.hh"

#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

// Hadrons
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4MesonConstructor.hh"

// Bosons
#include "G4BosonConstructor.hh"

// Leptons
#include "G4LeptonConstructor.hh"

// Other Particles
#include "G4ShortLivedConstructor.hh"

// Process options
#include "G4LossTableManager.hh"

// Physics List Helper
#include "G4PhysicsListHelper.hh"

#include "G4StepLimiter.hh"

#include <set>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSPhysicsList::TSPhysicsList()
{
  defaultCutValue = 1. * CLHEP::mm;

  fConstructors.push_back(fEmPhysics_opt4);
  fConstructors.push_back(fDecayPhysics);
  fConstructors.push_back(fRadDecayPhysics);
  fConstructors.push_back(fHadronInelasticPhysics);
  fConstructors.push_back(fHadronElasticPhysics);
  fConstructors.push_back(fIonElasticPhysics);
  fConstructors.push_back(fIonBinaryCascadePhysics);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSPhysicsList::~TSPhysicsList()
{
  for(auto ite : fConstructors)
    delete ite;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSPhysicsList::ConstructParticle()
{
  for(auto c : fConstructors)
  {
    c->ConstructParticle();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSPhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();

  for(auto c : fConstructors)
  {
    c->ConstructProcess();
  }

  std::set<G4String> step_limit_particles;
  // standard particles
  step_limit_particles.insert("e-");
  step_limit_particles.insert("e+");
  step_limit_particles.insert("alpha");
  step_limit_particles.insert("He3");
  step_limit_particles.insert("GenericIon");
  step_limit_particles.insert("proton");
  step_limit_particles.insert("neutron");
  // more ~exotic particles
  step_limit_particles.insert("pi+");
  step_limit_particles.insert("pi-");
  step_limit_particles.insert("mu+");
  step_limit_particles.insert("mu-");

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  while((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();
    G4String pname                 = particle->GetParticleName();

    if(step_limit_particles.find(pname) != step_limit_particles.end() ||
       particle->GetPDGCharge())
    {
      ph->RegisterProcess(new G4StepLimiter, particle);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
