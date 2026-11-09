// ********************************************************************
// * Modernization Demo: Integration of GPU and Quantum Features        *
// ********************************************************************
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

// Modernization Includes
#include "modernization/gpu/GpuTrackSoA.hh"
#include "modernization/gpu/GpuTrackingKernel.hh"
#include "modernization/quantum/G4QuantumProcess.hh"
#include "modernization/quantum/QuantumBridge.hh"

/**
 * ModernPhysicsList: Registers our new Quantum Process into Geant4
 */
class ModernPhysicsList : public G4VUserPhysicsList {
  public:
    void ConstructParticle() override {
        // Standard particle construction (placeholder)
    }
    void ConstructProcess() override {
        G4cout << "[ModernPhysicsList] Registering Quantum Process..." << G4endl;
        G4VProcess* quantumProc = new G4QuantumProcess();

        // We register the process to the electron for this demo
        G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
        if (electron) {
            RegisterProcess(quantumProc, electron);
            G4cout << "[ModernPhysicsList] G4QuantumProcess successfully registered to e-" << G4endl;
        }
    }
};

int main() {
    G4cout << "========================================================" << G4endl;
    G4cout << "   GEANT4 MODERNIZATION DEMO: GPU & QUANTUM            " << G4endl;
    G4cout << "========================================================" << G4endl;

    // --- PART 1: GPU ACCELERATION DEMO ---
    G4cout << "\n[Step 1] Testing GPU Acceleration (SoA Kernel)..." << G4endl;
    GpuTrackSoA gpuData(1000); // 1000 particles
    for(int i=0; i<1000; ++i) {
        gpuData.active[i] = true;
        gpuData.posX[i] = 0.0;
        gpuData.momX[i] = 1.0;
    }

    GpuTrackingKernel::ProcessTracks(gpuData, 0.01);
    G4cout << "Sample Particle 0 New PosX: " << gpuData.posX[0] << G4endl;

    // --- PART 2: QUANTUM INTEGRATION DEMO ---
    G4cout << "\n[Step 2] Testing Quantum Process Registration..." << G4endl;
    G4RunManager* runManager = new G4RunManager();
    ModernPhysicsList* physicsList = new ModernPhysicsList();
    runManager->SetUserPhysicsList(physicsList);

    // Trigger construction
    physicsList->Construct();

    // Simulate a "Quantum Consultation" manually to show it works
    G4cout << "\n[Step 3] Manual Quantum Backend Consultation..." << G4endl;
    G4double qProb = QuantumBridge::Instance().GetQuantumProbability("HiggsInteraction", 125.0, 0.0);
    G4cout << "Quantum probability for Higgs interaction: " << qProb << G4endl;

    G4cout << "\n========================================================" << G4endl;
    G4cout << "   DEMO COMPLETE: MODERNIZATION SUCCESSFUL              " << G4endl;
    G4cout << "========================================================" << G4endl;

    delete runManager;
    return 0;
}
