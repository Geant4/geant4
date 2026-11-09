// ********************************************************************
// * Quantum Integration Prototype: QuantumBridge                      *
// ********************************************************************
#ifndef QUANTUM_BRIDGE_HH
#define QUANTUM_BRIDGE_HH

#include "globals.hh"
#include <string>

/**
 * QuantumBridge provides the interface between Geant4 (C++) and
 * Quantum Computing Frameworks (Python/Qiskit/Cirq).
 */
class QuantumBridge {
  public:
    static QuantumBridge& Instance() {
        static QuantumBridge instance;
        return instance;
    }

    // Simulates a quantum circuit call to determine interaction probability
    G4double GetQuantumProbability(const std::string& processName, double energy, double momentum) {
        // Prototype: returning a dummy value.
        // Future: This will use pybind11 to call a Qiskit circuit.
        G4cout << "[QuantumBridge] Consulting quantum backend for " << processName << "..." << G4endl;
        return 0.015; // Example probability
    }

  private:
    QuantumBridge() {}
    QuantumBridge(const QuantumBridge&) = delete;
    void operator=(const QuantumBridge&) = delete;
};

#endif
