// ********************************************************************
// * Quantum Integration: G4QuantumProcess Implementation             *
// ********************************************************************
#include "G4QuantumProcess.hh"
#include "QuantumBridgePy.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

G4QuantumProcess::G4QuantumProcess() : G4VProcess("QuantumProcess", fNotDefined) {}

G4QuantumProcess::~G4QuantumProcess() {}

G4VParticleChange* G4QuantumProcess::AlongStepDoIt(const G4Track& track, const G4Step& stepData) {
    G4double prob = ConsultQuantumBackend(track);

    // If the quantum probability triggers an interaction
    if (G4UniformRand() < prob) {
        G4cout << "[QuantumProcess] Interaction triggered by Quantum Backend for Track ID: "
               << track.getTrackID() << G4endl;

        // For prototype, we return a basic particle change (e.g., energy loss)
        G4ParticleChange* change = new G4ParticleChange();
        change->SetEnergyChange(-0.1 * MeV); // Dummy energy loss
        return change;
    }

    return nullptr; // No interaction
}

G4VParticleChange* G4QuantumProcess::PostStepDoIt(const G4Track& track, const G4Step& stepData) {
    return nullptr;
}

G4VParticleChange* G4QuantumProcess::AtRestDoIt(const G4Track& track, const G4Step& stepData) {
    return nullptr;
}

G4double G4QuantumProcess::AlongStepGetPhysicalInteractionLength(const G4Track& track,
                                                                   G4double previousStepSize,
                                                                   G4double currentMinimumStep,
                                                                   G4double& proposedSafety,
                                                                   G4GPILSelection* selection) {
    // Return a dummy interaction length for prototype
    return 10.0 * cm;
}

G4double G4QuantumProcess::ConsultQuantumBackend(const G4Track& track) {
    return QuantumBridgePy::Instance().GetQuantumProbability(
        "QuantumScattering",
        track.GetKineticEnergy(),
        track.GetMomentum().mag()
    );
}
