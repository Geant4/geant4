\page ExampleChem1 Example chem1

## General description
  
  This example is provided by the Geant4-DNA collaboration.

  These processes and models are further described at:
  http://geant4-dna.org

  Any report or published results obtained using the Geant4-DNA software shall 
  cite the following Geant4-DNA collaboration publications:
  Phys. Med. 31 (2015) 861-874
  Med. Phys. 37 (2010) 4692-4708

  How to activate chemistry code.

## GEOMETRY DEFINITION
 
  It is a simple box which represents a 'semi infinite' homogeneous medium.
 	
  Two parameters define the geometry :
   - the material of the box -- for Geant4-DNA it has to be water.
   - the full size of the box.
	 	
  The default geometry is constructed in DetectorConstruction class.
 	
## PHYSICS LIST
  
  PhysicsList is Geant4 modular physics list using G4EmDNAPhysics & 
  G4EmDNAChemistry constructors.
 	 
## AN EVENT: THE PRIMARY GENERATOR
 
  The primary kinematic consists of a single particle starting at the center of 
  the box. The type of the particle and its energy are set in the 
  PrimaryGeneratorAction class, and can be changed via the G4 build-in commands 
  of G4ParticleGun class.
  The chemistry module is triggered in the StackingAction class when all 
  physical tracks have been processed.

## OUTPUT

  Physics initialization and the defined reaction table are printed.
  G4ITStepManager processes the chemical stage time step after time step.
  Chemical reactions are printed.

## HOW TO START

  To run the example in batch mode

      ./chem1 -mac beam.in

  or

      ./chem1

  then the macro beam.in is processed by default

  In interactive mode, run:

      ./chem1 -gui

  or

      ./chem1 -gui gui.mac
