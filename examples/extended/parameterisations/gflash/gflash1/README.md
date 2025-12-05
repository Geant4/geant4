\page Examplegflash1 Example gflash1

## Geometry and SD

This example uses only the mass geometry, with each crystal defined as a volume,
with parametrisation attached to the envelope in the mass geometry.
Geometry and sensitive detector are defined in:
- ExGflash1DetectorConstruction
- ExGflash1SensitiveDetector

## Details of implementation

To use GFLASH the user has to implement the following:

  - ExGflash1DetectorConstruction::ConstructSDandField() : \n
    Here GFLASH has to be initialized and assigend to the envelope,
    where it should be active (here our calorimeter = caloLog )
    ```cpp
    // **********************************************
    // * Initializing shower modell
    // ***********************************************
    G4cout << "Creating shower parameterization models" << G4endl;
    fFastShowerModel = new GFlashShowerModel("fFastShowerModel", fRegion);
    fParameterisation = new GFlashHomoShowerParameterisation(pbWO4);
    fFastShowerModel->SetParameterisation(*fParameterisation);
    fParticleBounds = new GFlashParticleBounds();
    fFastShowerModel->SetParticleBounds(*fParticleBounds);
    fHitMaker = new GFlashHitMaker();
    fFastShowerModel->SetHitMaker(*fHitMaker);
    G4cout<<"end shower parameterization."<<G4endl;
    // **********************************************
    ```

  - ExGflash1SensitiveDetector: \n
    It is mandatory to use G4VGFlashSensitiveDetector as (additional)
    base class for the sensitive detector.
    Here it is necessary to implement a seperate
    interface, where the GFlash spots are processed.
    (ProcessHits(G4GFlashSpot*aSpot ,G4TouchableHistory* ROhist))
    The separate interface is used, because the GFLASH spots contains
    (naturally) less information than the full simulation.

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
