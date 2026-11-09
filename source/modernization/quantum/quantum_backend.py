# quantum_backend.py
# This script is called by Geant4 via pybind11
try:
    from qiskit import QuantumCircuit, Aer, execute
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False

def calculate_interaction_probability(process_name, energy, momentum):
    """
    Simulates a quantum circuit to determine the probability
    of a particle interaction.
    """
    if not QISKIT_AVAILABLE:
        print("[Python] Qiskit not installed! Returning classical fallback.")
        return 0.01 # Fallback value

    # 1. Create a Quantum Circuit with 1 qubit
    qc = QuantumCircuit(1)

    # 2. Use the energy to determine the rotation angle (Quantum State)
    # We map energy (MeV) to a rotation between 0 and Pi
    import numpy as np
    theta = (energy % 100) * (np.pi / 100)
    qc.rx(theta, 0)

    # 3. Measure the qubit
    qc.measure_all()

    # 4. Execute the circuit on a local simulator
    backend = Aer.get_backend('qasm_simulator')
    job = execute(qc, backend, shots=1024)
    result = job.result()
    counts = result.get_counts()

    # Calculate probability of state |1> (the "Interaction" state)
    total_shots = sum(counts.values())
    interaction_count = counts.get('1', 0)

    return interaction_count / total_shots

print("[Python] Quantum Backend Initialized successfully.")
