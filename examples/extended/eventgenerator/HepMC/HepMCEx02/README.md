\page ExampleHepMCEx02 Example HepMCEx02

  This example demonstrates how to interface primary particles in Geant4
with various event generators via the HepMC Monte Carlo event interface.
This is another example having the same generator action as 
[HepMCEx01](../../html_HepMCEx01/html/ExampleHepMCEx01.html), 
but much simpler user control.

## Primary Generator

 H02PrimaryGeneratorAction has HepMCG4Interface as the generator.
There are two types of generators provided as samples. One generator reads 
primary information from a HepMC Ascii file (data/example_MyPythia.dat).
The other one generates primaries directly invoking PYTHIA routines 
in every event.

## Geometry

  A simplified collider-type geometry, which consists of 
    - endcap calorimeter (a set of tubes filled with lead), 
    - barrel calorimeter (tube filled with lead),
    - barrel muon detector (8 sets of plates filled with Ar),
    - endcap muon detecror, (a set of tubes filled with Ar) and
    - uniform magnetic field along the z axis of 3 Tesla at the 
      central region.

## Physics List

  FTFP_BERT predefined physics list

## User actions

  All particles except muons are killed in the calorimeter section.

## Installation

 See ["eventgenerator/HepMC" page](../../html/Examples_HepMC.html) how to build this example.
  
## Execution
```
% /HepMCEx02 hepmc_pygen.in
```
