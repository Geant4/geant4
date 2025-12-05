\page Examplemicroprox Example microprox

Author: S. Incerti et al. \n
Date: March 2nd, 2019     \n
Email: incerti@lp2ib.in2p3.fr

(c) The Geant4-DNA collaboration.

This example shows how to compute proximity functions
in liquid water using exclusively Geant4-DNA
physics processes and models.

This example is provided by the Geant4-DNA collaboration.

These processes and models are further described at:
http://geant4-dna.org

Any report or published results obtained using the Geant4-DNA software shall
cite the following Geant4-DNA collaboration publications: \n
J. Appl. Phys. (2019) in press  \n
Med. Phys. 51 (2024) 5873–5889  \n
Med. Phys. 45 (2018) e722-e739  \n
Phys. Med. 31  (2015) 861-874   \n
Med. Phys. 37  (2010) 4692-4708 \n
Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178

## Geometry

An infinite box of liquid water.

## Incident particles

Particles can be selected from the microprox.in macro
as well as their incident energy.
They are shot from the center of the box.
Tracking cut can also be selected (as energy).

## Physics

The default Geant4-DNA physics constructor 2 is used in
the PhysicsList class. Alternative constructor can be
selected from microprox.in

## Scoring of enery deposition

Energy depositions are scored in spherical shells from randomly selected hits.
The user can select the dimensions of the shells as well as radius steps in TrackerSD.

## Run

The code can be run using:
```
./microprox microprox.in
```

## Results

Results can be analyzed after the run using:
```
root plot.C
```

The distribution of t is shown by default.
