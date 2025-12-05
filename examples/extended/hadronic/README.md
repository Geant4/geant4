\page Examples_hadronic Category "hadronic"

 Examples in this directory demonstrate specific hadronic physics simulation 
 with histogramming.

\ref ExampleHadr00

This example demonstrates a usage of G4PhysListFactory to build 
Physics List and G4HadronicProcessStore to access cross sections.

\ref ExampleHadr01

This example application is based on the application IION developed for
simulation of proton or ion beam interaction with a water target. Different 
aspects of beam target interaction are demonstrating in the example including 
longitudinal profile of energy deposition, spectra of secondary  particles,
spectra of particles leaving the target. 

\ref ExampleHadr02

This example application is providing simulation of ion beam interaction with different 
targets. Hadronic aspects of beam target interaction are demonstrated in the example 
including longitudinal profile of energy deposition, spectra of secondary  particles,
isotope production spectra. 

\ref ExampleHadr03

This example demonstrates how to compute total cross section from the direct evaluation of the 
mean free path ( see below, item Physics), how to identify nuclear reactions, how to plot 
energy spectrum of secondary particles.	

\ref ExampleHadr04

This example is focused on neutronHP physics, especially neutron transport,
including thermal scattering.
See A.R. Garcia, E. Mendoza, D. Cano-Ott presentation at G4 Hadronic group 
meeting (04/2013) and note on G4NeutronHP package

\ref ExampleHadr05

Examples of hadronic calorimeters

\ref ExampleHadr06

This example demonstrates survey of energy deposition and particle's flux from 
a hadronic cascade.

\ref ExampleHadr07

Survey energy deposition and particle's flux from an hadronic cascade.
Use PhysicsConstructor objects rather than predefined G4 PhysicsLists.
Show how to plot a depth dose profile in a rectangular box.

\ref ExampleHadr08

This example shows how to get "hadronic model per region" using generic
biasing: in particular, it is shown how to use "FTFP+INCLXX" in one region,
while using the default "FTFP+BERT" in all other regions. 
Notice that we use the generic biasing machinery, but the actual weights
of all tracks remain to the usual value (1.0) as in the normal (unbiased)
case.

\ref ExampleHadr09

This example shows how to use Geant4 as a generator for simulating
inelastic hadron-nuclear interactions.
Notice that the Geant4 run-manager is not used.

\ref ExampleHadr10

This example aims to test the treatment of decays in Geant4.
In particular, we want to test the decays of the tau lepton, charmed and
bottom hadrons, and the use of pre-assigned decays.

\ref ExampleFissionFragment

This example demonstrates the Fission Fragment model as used within the
neutron_hp model. It will demonstrate the capability for fission product
containment by the cladding in a water moderated sub-critical assembly. It could
also be further extended to calculate the effective multiplication factor of
the subcritical assembly for various loading schemes.

\ref Examples_FlukaCern

A set of 2 examples, demonstrating how to make use of 
the interface to `FLUKA` hadron-nucleus inelastic physics in a G4 application. \n
The examples are at the process (interaction) level, but a physics list 
(G4_HP_CernFLUKAHadronInelastic_PhysicsList) is also made available. <br>
The interface to `FLUKA` itself is also included.

\ref ExampleNeutronSource

NeutronSource is an example of neutrons production. It illustrates the cooperative work
of nuclear reactions and radioactive decay processes.
It survey energy deposition and particle's flux.
It uses PhysicsConstructor objects.

\ref Examples_ParticleFluence

This example aims to monitor the particle fluence for various particle types
and set-ups. The particle fluence at a given position is defined as the
average number of particles crossing a unit surface in that position
(normalized per one incident primary). The particle fluence is conveniently
estimated by summing the particles' track lengths in a thin scoring volume
and dividing for the cubic volume of such a scoring volume.

