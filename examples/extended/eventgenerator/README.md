\page Examples_eventgenerator Category "eventgenerator"

 Examples in this directory demonstrate various ways of primary event
 generation.

\ref ExampleparticleGun

This example demonstrates 4 ways of the usage of G4ParticleGun shooting
primary particles in different cases.

\ref Exampleexgps

This example demonstrates the usage of G4GeneralParticleSource for generating 
primary incident particle according to user defined distributions.

\ref ExampleuserPrimaryGenerator

This example shows how to create a primary event including several vertices and 
several primary particles per vertex.

\ref Examples_HepMC

This directory contains examples for using HepMC as an interface with 
various Monte Carlo event generators, such as PYTHIA.
It also include an example for demonstrating MC truth handling with HepMC. 

\ref Examples_pythia

This directory contains the following examples:

a) use of Pythia6 as Monte Carlo event generator, interfaced with Geant4, 
and showing how to implement an external decayer based on Pythia6.
The feature is activated by setting environment variable PYTHIA6 to point
to the Pythia6 installation area.
For details, please see \ref Exampledecayer6.

b) use of Pythia8 as an external decayer to replace native Geant4 decay
tables for such resonances as tau+/- and B+/-, and to supplement Pythia8-based
decay tables to those resonances where Geant4 native decay features are not
implemented. 
For details, please see \ref Examplepy8decayer.
