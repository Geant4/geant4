\page Examplegflash3 Example gflash3

## Geometry and SD

This example uses mass geometry to create the main volume (homogeneous material) and use it
as an envelope for the parametrisation, but the readout geometry (crystals)
are defined in the parallel geometry, together with the sensitive detector
to store the hits.
Geometry and sensitive detector are defined in:
- ExGflash3DetectorConstruction
- ExGflash3ParallelWorld
- ExGflash3SensitiveDetector

## Details of implementation:

- ExGflash3.cc:
   Parallel world needs to be registered;
   Physics of the parallel world needs to be registered so sensitive detector can
   collect hits;
   Fast simulation is activated for mass world (where envelope is);

- ExGflash3DetectorConstruction:
   Only main volume (box) with material is created;
   Creation of G4Region associated to G4LogicalVolume of that box;
   Initialization of GFlash, attaching it to the envelope (G4Region);

- ExGflash3ParallelWorld:
   Construction of identical volume for the main box as in the mass geometry,
   but with dummy material (it is not used anyway);
   Construction of individual crystals for the readout geometry;
   Creation of the sensitive detector;

- ExGflash3SensitiveDetector:
   Uses pointer to ExGflash3DetectorConstruction to get the crystals for the readout;

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

See also [Category "parameterisations/gflash"](../../html/Examples_gflash.html) documentation.

