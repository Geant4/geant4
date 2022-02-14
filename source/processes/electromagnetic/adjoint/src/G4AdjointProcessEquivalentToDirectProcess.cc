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
// ------------------------------------------------------------
//
// This class is for adjoint process equivalent to direct process
//
// ------------------------------------------------------------
//   Created by L.Desorgher          25 Sept. 2009
// ------------------------------------------------------------

#include "G4AdjointProcessEquivalentToDirectProcess.hh"

#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

G4AdjointProcessEquivalentToDirectProcess::
  G4AdjointProcessEquivalentToDirectProcess(
    const G4String& aName, G4VProcess* aProcess,
    G4ParticleDefinition* fwd_particle_def)
  : G4VProcess(aName)
{
  fDirectProcess  = aProcess;
  theProcessType  = fDirectProcess->GetProcessType();
  fFwdParticleDef = fwd_particle_def;
}

G4AdjointProcessEquivalentToDirectProcess::
  ~G4AdjointProcessEquivalentToDirectProcess()
{
  if(fDirectProcess != nullptr)
    delete fDirectProcess;
}

void G4AdjointProcessEquivalentToDirectProcess::
  ResetNumberOfInteractionLengthLeft()
{
  fDirectProcess->ResetNumberOfInteractionLengthLeft();
}

G4double G4AdjointProcessEquivalentToDirectProcess::
  AlongStepGetPhysicalInteractionLength(const G4Track& track,
                                        G4double previousStepSize,
                                        G4double currentMinimumStep,
                                        G4double& proposedSafety,
                                        G4GPILSelection* selection)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  // Call the direct process
  G4double GPIL = fDirectProcess->AlongStepGetPhysicalInteractionLength(
    track, previousStepSize, currentMinimumStep, proposedSafety, selection);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return GPIL;
}

G4double
G4AdjointProcessEquivalentToDirectProcess::AtRestGetPhysicalInteractionLength(
  const G4Track& track, G4ForceCondition* condition)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  // Call the direct process
  G4double GPIL =
    fDirectProcess->AtRestGetPhysicalInteractionLength(track, condition);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return GPIL;
}

G4double
G4AdjointProcessEquivalentToDirectProcess::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double previousStepSize, G4ForceCondition* condition)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());

  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  // Call the direct process
  G4double GPIL = fDirectProcess->PostStepGetPhysicalInteractionLength(
    track, previousStepSize, condition);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return GPIL;
}

G4VParticleChange* G4AdjointProcessEquivalentToDirectProcess::PostStepDoIt(
  const G4Track& track, const G4Step& stepData)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());

  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  // Call the direct process
  G4VParticleChange* partChange = fDirectProcess->PostStepDoIt(track, stepData);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return partChange;
}

G4VParticleChange* G4AdjointProcessEquivalentToDirectProcess::AlongStepDoIt(
  const G4Track& track, const G4Step& stepData)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());

  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  // Call the direct process
  G4VParticleChange* partChange =
    fDirectProcess->AlongStepDoIt(track, stepData);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return partChange;
}

G4VParticleChange* G4AdjointProcessEquivalentToDirectProcess::AtRestDoIt(
  const G4Track& track, const G4Step& stepData)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());

  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  // Call the direct process
  G4VParticleChange* partChange = fDirectProcess->AtRestDoIt(track, stepData);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return partChange;
}

G4bool G4AdjointProcessEquivalentToDirectProcess::IsApplicable(
  const G4ParticleDefinition&)
{
  return fDirectProcess->IsApplicable(*fFwdParticleDef);
}

void G4AdjointProcessEquivalentToDirectProcess::BuildPhysicsTable(
  const G4ParticleDefinition&)
{
  return fDirectProcess->BuildPhysicsTable(*fFwdParticleDef);
}

void G4AdjointProcessEquivalentToDirectProcess::PreparePhysicsTable(
  const G4ParticleDefinition&)
{
  return fDirectProcess->PreparePhysicsTable(*fFwdParticleDef);
}

G4bool G4AdjointProcessEquivalentToDirectProcess::StorePhysicsTable(
  const G4ParticleDefinition*, const G4String& directory, G4bool ascii)
{
  return fDirectProcess->StorePhysicsTable(fFwdParticleDef, directory, ascii);
}

G4bool G4AdjointProcessEquivalentToDirectProcess::RetrievePhysicsTable(
  const G4ParticleDefinition*, const G4String& directory, G4bool ascii)
{
  return fDirectProcess->RetrievePhysicsTable(fFwdParticleDef, directory,
                                              ascii);
}

void G4AdjointProcessEquivalentToDirectProcess::StartTracking(G4Track* track)
{
  // Change the particle definition to the direct one
  G4DynamicParticle* theDynPart =
    const_cast<G4DynamicParticle*>(track->GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();

  G4DecayProducts* decayProducts =
    const_cast<G4DecayProducts*>(theDynPart->GetPreAssignedDecayProducts());
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*) (0));
  theDynPart->SetDefinition(fFwdParticleDef);

  fDirectProcess->StartTracking(track);

  // Restore the adjoint particle definition to the direct one
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);

  return;
}

void G4AdjointProcessEquivalentToDirectProcess::EndTracking()
{
  fDirectProcess->EndTracking();
}
