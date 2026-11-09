// ********************************************************************
// * Quantum Integration: QuantumBridgePy Wrapper                     *
// ********************************************************************
#ifndef QUANTUM_BRIDGE_PY_HH
#define QUANTUM_BRIDGE_PY_HH

#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include "globals.hh"
#include <string>

namespace py = pybind11;

/**
 * QuantumBridgePy handles the embedding of the Python interpreter
 * within the Geant4 process to communicate with Qiskit.
 */
class QuantumBridgePy {
  public:
    static QuantumBridgePy& Instance() {
        static QuantumBridgePy instance;
        return instance;
    }

    // Initialize the Python interpreter and load the backend script
    void Initialize() {
        try {
            if (!py::module_exists("qiskit")) {
                G4cerr << "[QuantumBridgePy] Warning: Qiskit not found in Python environment!" << G4endl;
            }
            // Start the interpreter
            py::initialize_interpreter();
            // Import the backend script we created
            py::module_ backend = py::module_::import("quantum_backend");
            G4cout << "[QuantumBridgePy] Python Quantum Backend loaded successfully." << G4endl;
        } catch (const py::error_already_set& e) {
            G4cerr << "[QuantumBridgePy] Python Error: " << e.what() << G4endl;
        }
    }

    // Call the Python function: calculate_interaction_probability
    G4double GetQuantumProbability(const std::string& processName, double energy, double momentum) {
        try {
            auto backend = py::module_::import("quantum_backend");
            auto result = backend.attr("calculate_interaction_probability")(processName, energy, momentum);
            return result.cast<double>();
        } catch (const py::error_already_set& e) {
            G4cerr << "[QuantumBridgePy] Runtime Python Error: " << e.what() << G4endl;
            return 0.0; // Fail-safe
        }
    }

    void Finalize() {
        py::finalize_interpreter();
    }

  private:
    QuantumBridgePy() {}
    QuantumBridgePy(const QuantumBridgePy&) = delete;
    void operator=(const QuantumBridgePy&) = delete;
};

#endif
