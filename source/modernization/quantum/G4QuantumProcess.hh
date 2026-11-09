// ********************************************************************
// * Quantum Integration Prototype: G4QuantumProcess                   *
// ********************************************************************
#ifndef G4_QUANTUM_PROCESS_HH
#define G4_QUANTUM_PROCESS_HH

#include "G4VProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"

/**
 * G4QuantumProcess is a physics process that delegates its calculations
 * to a Quantum Backend instead of using classical cross-sections.
 */
class G4QuantumProcess : public G4VProcess {
  public:
    G4QuantumProcess();
    virtual ~G4QuantumProcess();

    // Core Physics Overrides
    virtual G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& stepData) override;
    virtual G4VParticleChange* AlongStepDoIt(const G4Track& track, const G4Step& stepData) override;
    virtual G4VParticleChange* AtRestDoIt(const G4Track& track, const G4Step& stepData) override;

    virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track& track,
                                                           G4double previousStepSize,
                                                           G4double currentMinimumStep,
                                                           G4double& proposedSafety,
                                                           G4GPILSelection* selection) override;

  private:
    // Internal method to call the Quantum Bridge
    G4double ConsultQuantumBackend(const G4Track& track);
};

#endif
