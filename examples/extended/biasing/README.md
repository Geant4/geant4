\page Examples_biasing Category "biasing"

## B01, B02 and B03

B01, B02 and B03 applications demonstrate the usage of different variance
reduction techniques supported in Geant4, or possible from the user
applications.

### General remark to variance reduction

The tools provided for importance sampling (or geometrical splitting and
Russian roulette) and for the weight window technique require the user to 
have a good understanding of the physics in the problem. This is because 
the user has to decide which particle types have to be biased, define the 
cells (physical volumes, replicas) and assign importances or weight 
windows to that cells. If this is not done properly it can not be 
expected that the results describe a real experiment. The examples given 
here only demonstrate how to use the tools technically. They don't intend 
to produce physical correct results.

### General remark to scoring

Scoring is carried out using the built-in Multifunctional detectors. For
parallel geometries this requires a special scoring physics process.
See examples/extended/runAndEvent (especailly RE05) for clarification.

### Known problems - should not happen

In the following scenario it can happen that a particle is not
biased and it's weight is therefore not changed even if it crosses
a boundary where biasing should happen.
Importance and weight window sampling create particles on boundaries 
between volumes. If the GPIL method of a physical process returns 
0 as step length for a particle on a boundary and if the PostStepDoIt of
that process changes the direction of the particle to go back in the 
former volume the biasing won't be invoked. 
This will produce particles with weights that do not correspondent to the
importance of the current volumes.

### Further information:

Short description of importance sampling and scoring:
https://geant4.web.cern.ch/collaboration/working_groups/geometryTransport/#development-documents (Under the Event Biasing & Tallies Section)

### \ref ExampleB01

The example uses importance sampling or the weight window technique 
according to an input parameter.

### \ref ExampleB02

This example uses a parallel geometry to define G4GeometryCell objects
for scoring and importance sampling. The output should be equivalent to B01.

A modular approach is applied to the physicslist and the extension for biasing.
The parallel geometry is included in this extension.

### \ref ExampleB03

This example uses a parallel geometry to define G4GeometryCell objects
for scoring and importance sampling. The output should be statistically 
equivalent to B02 (and B01).

This demonstrates a customised "flat" physics implementation with the addition
of biasing. Complementary approach to the modular physics lists of B01 and B02


## Generic biasing examples GB01 - GB06

These examples illustrate the usage of a biasing scheme implemented since
version Geant4 10.0.
The scheme is meant to be extensible, not limited to these six examples.

### \ref ExampleGB01

This example illustrates how to bias process cross-sections in this scheme.

### \ref ExampleGB02

Illustrates a force collision scheme similar to the MCNP one.

### \ref ExampleGB03

Illustrates geometry based biasing.

### \ref ExampleGB04

Illustrates a bremsstrahlung splitting.

### \ref ExampleGB05

Illustrates a "splitting by cross-section" technique: a splitting-based
technique using absorption cross-section to control the neutron population.

### \ref ExampleGB06

Illustrates the usage of parallel geometries with generic biasing.

### \ref ExampleGB07

Illustrates the usage of leading particle biasing with generic biasing.

## Reverse MonteCarlo Technique example

### \ref ExampleReverseMC01

Example illustrating the use of the Reverse Monte Carlo (RMC) mode in a Geant4
application. See details in \ref ExampleReverseMC01 Example README page
\endlink.
