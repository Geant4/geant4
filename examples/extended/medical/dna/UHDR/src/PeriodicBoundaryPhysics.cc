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

#include "PeriodicBoundaryPhysics.hh"

#include "PeriodicBoundaryProcess.hh"

#include "G4PhysicsConstructorFactory.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4_DECLARE_PHYSCONSTR_FACTORY(PeriodicBoundaryPhysics);
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PeriodicBoundaryPhysics::PeriodicBoundaryPhysics(const G4String& name, G4bool per_x, G4bool per_y,
                                                 G4bool per_z)
  : G4VPhysicsConstructor(name), fPeriodicX(per_x), fPeriodicY(per_y), fPeriodicZ(per_z)
{
  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PeriodicBoundaryPhysics::ConstructParticle() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PeriodicBoundaryPhysics::ConstructProcess()
{
  if (verboseLevel > 0) G4cout << "Constructing cyclic boundary physics process" << G4endl;

  auto* pbc =
    new PeriodicBoundaryProcess("Cyclic", fNotDefined, fPeriodicX, fPeriodicY, fPeriodicZ);

  if (verboseLevel > 0) {
    pbc->SetVerboseLevel(verboseLevel);
  }

  auto aParticleIterator = GetParticleIterator();

  aParticleIterator->reset();

  G4ProcessManager* processManager = nullptr;

  while ((*aParticleIterator)()) {
    G4ParticleDefinition* particle = aParticleIterator->value();

    G4String particleName = particle->GetParticleName();

    processManager = particle->GetProcessManager();

    if (!processManager) {
      ThrowException(particleName);
      return;
    }

    AddDiscreteProcess(pbc, *particle, processManager);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PeriodicBoundaryPhysics::ThrowException(const G4String& particleName)
{
  std::ostringstream o;
  o << "Particle " << particleName << "without a Process Manager";
  G4Exception("G4PeriodicBoundaryPhysics::ConstructProcess()", "", FatalException, o.str().c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PeriodicBoundaryPhysics::AddDiscreteProcess(PeriodicBoundaryProcess* periodicBoundaryProcess,
                                                 G4ParticleDefinition& particle,
                                                 G4ProcessManager* processManager)
{
  if (periodicBoundaryProcess->IsApplicable(particle)) {
    processManager->AddDiscreteProcess(periodicBoundaryProcess);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
