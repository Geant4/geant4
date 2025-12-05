\page Examples_gflash Category "parameterisations/gflash"

These examples demonstrate the use of the GFLASH parameterisation library.
They use the GFLASH equations (hep-ex/0001020, Grindhammer & Peters)
to parametrise electromagnetic showers in matter.
In these examples the calorimeter is a simple cube,
which consists of 10 x 10 crystals of PbWO4 (CMS like).

Briefly, whenever a e-/e+  enters the calorimeter, it is parametrised if it
has a minimum energy and the shower is expected to be contained
in the calorimeter (so called " parameterisation envelope").
If this is fullfilled the particle is killed, as well as all secondaries,
and the energy is deposited according to the GFLASH equations.

The examples show how to interface GFLASH to your application.
The simulation time is measured, so the user can see immediately
the speed up by using GFLASH.

Geometry and parametrisation is defined in different ways in the set of three equivalent
(in terms of produced showers) examples: \ref Examplegflash1, \ref Examplegflash2 and 
\ref Examplegflash3, to demonstrate
how to use the parametrisation, sensitive detectors and parallel geometry.
The classes which are the same in all three examples have the names with ExGflash prefix while
the names of classes specific to each example have the prexix ExGflash[1,2,3].

\ref Examplegflasha - allow histogramming of show profiles and fine tuning
of gflash parametrization for homogeneous medium.

Note: Instead of particle gun the gps class is used here for particle generation.

## Briefly

 Table below presents in which world/geometry (mass or parallel) each of the elements is defined.


|            Example           |  gflash1 |    gflash2   |    gflash3   |
|------------------------------|----------|--------------|--------------|
|       Block of material      | mass geo |   mass geo   |   mass geo   |
|  Crystals (readout geometry) | mass geo |   mass geo   | parallel geo |
|      Sensitive detector      | mass geo |   mass geo   | parallel geo |
| Envelope for parametrisation | mass geo | parallel geo |   mass geo   |


### \ref Examplegflash1

Uses only the mass geometry, with each crystal defined as a volume,
with parametrisation attached to the envelope in the mass geometry.
Geometry and sensitive detector are defined in:
- ExGflash1DetectorConstruction
- ExGflash1SensitiveDetector

### \ref Examplegflash2

Uses mass geometry to create volumes and to create a sensitive detector
for storing hits, but parametrisation is attached to the envelope
in the parallel geometry (see also examples/extended/parametrisations/Par01).
Geometry and sensitive detector are defined in:
- ExGflash2DetectorConstruction
- ExGflash2ParallelWorld
- ExGflash2SensitiveDetector

### \ref Examplegflash3

Uses mass geometry to create the main volume (homogeneous material) and use it
as an envelope for the parametrisation, but the readout geometry (crystals)
are defined in the parallel geometry, together with the sensitive detector
to store the hits.
Geometry and sensitive detector are defined in:
- ExGflash3DetectorConstruction
- ExGflash3ParallelWorld
- ExGflash3SensitiveDetector

## Macros

- vis.mac - macro for use in interactive mode (default, if no arguments are specified)
- test.mac - macro for tests: 50 GeV electrons are shot in the direction of the detector
           (along z axis), 10 times. As they enter the parametrisation envelope,
           the GFlash parametrisation is invoked and energy is deposited.
           The results are printed out:
            - energy in the most central crystal
            - energy in 3x3 crystals
            - energy in 5x5 crystals
            - number of created deposits
            - simulation time per event
