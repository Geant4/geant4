\page Examplegflash2 Example gflash2

## Geometry and SD

This example uses mass geometry to create the main volume (homogeneous material) and use it
as an envelope for the parametrisation, but the readout geometry (crystals)
are defined in the parallel geometry, together with the sensitive detector
to store the hits.
Geometry and sensitive detector are defined in:
- ExGflash2DetectorConstruction
- ExGflash2ParallelWorld
- ExGflash2SensitiveDetector

## Details of implementation

- ExGflash2.cc: \n
  Parallel world needs to be registered;
  Fast simulation is activated for parallel world (where envelope is);

- ExGflash2DetectorConstruction: \n
  Only main geometry and SD are created;

- ExGflash2ParallelWorld: \n
  Construction of identical volume for the main box as in the mass geometry,
  but with dummy material (it is not used anyway);
  Creation of G4Region associated to G4LogicalVolume;
  Initialization of GFlash, attaching it to the envelope (G4Region);

- ExGflash2SensitiveDetector: \n
  Uses pointer to ExGflash2ParallelWorld to get the crystals for the readout;

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

