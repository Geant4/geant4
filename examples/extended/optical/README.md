\page Examples_optical Category "optical"

This directory includes examples demonstrating the use of optical processes
in the simulation.

\ref ExampleOpNovice

Simulation of optical photons generation and transport.
Defines optical surfaces and exercises optical physics processes
(Cerenkov, Scintillation, Absorption, Rayleigh, ...). Uses stacking
mechanism to count the secondary particles generated.
Via the command line one can select an option to define the detector via a
gdml file. An example gdml file is provided that corresponds
to the  detector configuration defined in OpNoviceDetectorConstruction.cc.

\ref ExampleOpNovice2

Investigate optical properties and parameters. Details of optical
photon boundary interactions on a surface are recorded. Details 
of optical photon generation and transport are recorded.

\ref ExampleLXe

Multi-purpose detector setup implementing:
-# scintillation inside a bulk scintillator with PMTs
-# large wall of small PMTs opposite a Cerenkov slab to show the cone
-# plastic scintillator with wave-length-shifting fiber readout.

\ref Examplewls

This application simulates the propagation of photons inside a Wave Length
Shifting (WLS) fiber.
