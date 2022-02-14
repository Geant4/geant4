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
/////////////////////////////////////////////////////////////////////////////////
//  Class:    G4AdjointProcessEquivalentToDirectProcess
//  Author:         L. Desorgher
//  Organisation:   SpaceIT GmbH
//
//  Adjoint process equivalent to direct process, used for some multiple
//  scattering.
//  A virtual class for wrapper process objects.
/////////////////////////////////////////////////////////////////////////////////

#ifndef G4AdjointProcessEquivalentToDirectProcess_h
#define G4AdjointProcessEquivalentToDirectProcess_h 1

#include "G4VProcess.hh"

class G4AdjointProcessEquivalentToDirectProcess : public G4VProcess
{
 public:
  explicit G4AdjointProcessEquivalentToDirectProcess(
    const G4String& aName, G4VProcess* aProcess,
    G4ParticleDefinition* fwd_particle_def);

  ~G4AdjointProcessEquivalentToDirectProcess() override;

  G4VParticleChange* PostStepDoIt(const G4Track& track,
                                  const G4Step& stepData) override;

  G4VParticleChange* AlongStepDoIt(const G4Track& track,
                                   const G4Step& stepData) override;
  G4VParticleChange* AtRestDoIt(const G4Track& track,
                                const G4Step& stepData) override;

  G4double AlongStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4double currentMinimumStep, G4double& proposedSafety,
    G4GPILSelection* selection) override;

  G4double AtRestGetPhysicalInteractionLength(
    const G4Track& track, G4ForceCondition* condition) override;

  G4double PostStepGetPhysicalInteractionLength(
    const G4Track& track, G4double previousStepSize,
    G4ForceCondition* condition) override;

  G4bool IsApplicable(const G4ParticleDefinition&) override;

  void BuildPhysicsTable(const G4ParticleDefinition&) override;

  void PreparePhysicsTable(const G4ParticleDefinition&) override;

  G4bool StorePhysicsTable(const G4ParticleDefinition*,
                           const G4String& directory,
                           G4bool ascii = false) override;

  G4bool RetrievePhysicsTable(const G4ParticleDefinition*,
                              const G4String& directory,
                              G4bool ascii = false) override;

  void StartTracking(G4Track*) override;
  void EndTracking() override;

  void ResetNumberOfInteractionLengthLeft() override;

  G4AdjointProcessEquivalentToDirectProcess(G4AdjointProcessEquivalentToDirectProcess&) =
    delete;
  G4AdjointProcessEquivalentToDirectProcess& operator=(
    const G4AdjointProcessEquivalentToDirectProcess& right) = delete;

 private:
  G4ParticleDefinition* fFwdParticleDef;
  G4VProcess* fDirectProcess;
};

#endif
