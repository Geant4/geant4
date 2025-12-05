\page Exampledmparticle Example dmparticle

Geant4 application for testing of dark matter particles and processes
(based on TestEm8 and monopole examples)

This is very preliminary and simplified Geant4 example for light dark matter (LDM) particles. 
It consists of LDM bremsstrahlung process (in G4LDMPhysics constructor) 
for protons creating G4LDMPhotons. 
The latter decays creating G4LDMHi and G4LDMHiBar LDM scalar particles. They 
can scatter on nucleons and electrons. 

The first version has electromagnetic processes only. Hadron constructors can be added similarly 
to HadrN examples.

More LDM processes will be added according to experiment requirements

	
## GEOMETRY DEFINITION
 
 The geometry consists of a single block of a homogenous material,
 placed in a world.

 Parameters define the geometry :
 - the material of the box 
 - the thickness of the box 
 - the tranverse dimension of the box 

 The default is 130 cm of tungsten.
 Equivalent UI commands are following:
```
/testex/det/setMat G4_W
/testex/det/sizeX  130 cm
/testex/det/sizeYZ  30 cm
```

 The default geometry is constructed in DetectorConstruction class,
 but all of the above parameters can be changed interactively via
 the commands defined in the DetectorMessenger class. 

## PHYSICS LIST
 

 Physics Lists include EM standard physics and decay physics, additional 
 dark matter particle physics imlemented inside PhysicsList method. By
 default DMLPhoton is defined with the mass 0.5 GeV.
 To define different mass of this photon or enable DMLHi, DMLHiBar
 command line commands should be applied.

## THE PRIMARY GENERATOR

 The primary kinematic consists of a single particle which hits the
 block perpendicular to the input face. The type of the particle
 and its energy are set in the PrimaryGeneratorAction class, and can
 changed via the G4 build-in commands of ParticleGun class (see
 the macros provided with this example).
 The default is proton 100 GeV

## HOW TO START ?
 

 Execute Test  in 'batch' mode from macro files

```
% ./dmparticle dmparticle.in
% ./dmparticle dmparticle.in 0.2 0.4
```

 two extra numbers are masses in GeV of DMLPhoton and DMLHi

## HISTOGRAMS

 The result histogram:
 - energy deposition along the target

 The histogram is saved in Root file.
