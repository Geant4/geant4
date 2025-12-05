\page Examplemicroyz Example microyz

Author: S. Incerti et al.
Date: 17 Apr. 2015
Email: incerti@cenbg.in2p3.fr

(c) The Geant4-DNA collaboration.

This examples shows how to compute microdosimetry 
distributions y, z in liquid water using exclusively Geant4-DNA
physics processes and models.

This example is provided by the Geant4-DNA collaboration.

These processes and models are further described at:
http://geant4-dna.org

Any report or published results obtained using the Geant4-DNA software shall 
cite the following Geant4-DNA collaboration publications: \n
Phys. Med. 31 (2015) 861-874   \n
Med. Phys. 37 (2010) 4692-4708 \n
J. Appl. Phys. 122 (2017) 024303

## Geometry

A box of liquid water.

## Incident particles

Particles can be selected from the microyz.in macro
as well as their incident energy.
They are shot from the center of the box.

## Physics

The default Geant4-DNA physics constructor 2 is used in
the PhysicsList class. Alternative constructor can be
selected from the macro file microyz.in

Livermore and Penelope physics lists can be used as well.

Tracking cut and maximum step size can be selected in microyz.in

## Scoring of enery deposition

Energy depositions are scored in spheres randomly
placed along the incident particle track, 
using a weighted sampling.

The user can select the radius of the sphere in microyz.in using the command:
```
/microyz/det/Radius X unit
```
where X is the radius value and unit is specified.

## Run

The code can be run using:
```
./microyz microyz.in
```

## Results

Results can be analyzed after the run using:
```
root plot.C
```

The distribution of y is shown by default.

The following quantities are calculated: yF and yD.


