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
/// file: ParallelWorldPhysics.cc
/// brief: G4ParallelWorldPhysics, but without application to molecular species

#include "ParallelWorldPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(ParallelWorldPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParallelWorldPhysics::ParallelWorldPhysics(const G4String& name,
                                           G4bool layeredMass)
  : G4VPhysicsConstructor(name)
  , fLayeredMass(layeredMass)
{
  // ctor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParallelWorldPhysics::~ParallelWorldPhysics() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParallelWorldPhysics::ConstructParticle()
{
  // stub
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParallelWorldPhysics::ConstructProcess()
{
  // Make sure the parallel world registered
  G4TransportationManager::GetTransportationManager()->GetParallelWorld(
    namePhysics);

  // Add parallel world process//why do we need this
  auto theParallelWorldProcess = new G4ParallelWorldProcess(namePhysics);
  theParallelWorldProcess->SetParallelWorld(namePhysics);
  theParallelWorldProcess->SetLayeredMaterialFlag(fLayeredMass);

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();

  while((*particleIterator)())
  {
    G4ParticleDefinition* particle = particleIterator->value();
    if(!particle->IsShortLived() && particle->GetParticleType() != "Molecule"
       && theParallelWorldProcess->IsApplicable(*particle))
    {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      pmanager->AddProcess(theParallelWorldProcess);
      if(theParallelWorldProcess->IsAtRestRequired(particle))
      {
        pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9900);
      }
      pmanager->SetProcessOrderingToSecond(theParallelWorldProcess,
                                           idxAlongStep);
      pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9900);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
