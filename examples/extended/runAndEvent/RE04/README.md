\page ExampleRE04 Example RE04


 Contact : M.Asai (SLAC)

## Introduction

 This example demonstrates how to define a layered mass
geometry in parallel world. In the mass (tracking) world,
there are two boxes only. One is the world volume and the 
other is a box in the world. They both are made of air.
Thus, if tracks do not see materials (water and lead)
defined in the parallel world, they rarely interact.
In the parallel world, there are boxes made of water and
lead.

### Geometry

 RE04DetectorConstruction defines the mass (tracking)
geometry. It firstly defines all materials which apear
either in mass world or parallel world. Then in SetupGeometry()
method, it defines the world volume and a box named "phantom".
Both boxes are made of air.

 RE04ParallelWorldConstruction defines the parallel world.
For a parallel world, solid, logical and physical volumes
which represent parallel world must not be created here but
should be taken through G4VUserParallelWorld::GetWorld()
method which creates clones of solid, logical and physical
volumes of the world volume of the mass world. Please note
that this cloned logical volume of the parallel world volume
does not have a valid pointer to aa material but null.

 In the parallel world, if a logical volume has a valid
material pointer, a track in this volume (precisely saying
a physical volume which is made of this logical volume)
will see the material defined in this logical volume,
regardless of the material in the mass geometry. If a
logical volume has a null material pointer, a track will
see the ordinary material defined in the mass world.

 RE04ParallelWorldConstruction defines one placement
volume of box-shape, which is made of water, and a mother
box (placement volume with null material pointer), which
contains parameterized volumes. RE04ParallelWorldParam
class defines a parameterization of the parameterized
volume "paraPara", which represents two boxes at different
locations and made of water and lead respectively.

### Physics

 The physics list is taken from referenced physics-list FTFP_BERT
in Geant4.

## Macro files

 The macro file score.mac defines a scoring mesh which covers
the "Phantom" and scores energy deposition. It shoots 1000
primary particles (by default 10 GeV muon-). Though the mass
world has only air, given tracks, both primary muons and 
secondary particles see water and lead defined in the parallel 
world, you will see the energy deposition is not evenly
distributed.

The macro file batch.mac defines the same setup as score.mac,
but visualization is disabled and the score dump is activated.

## User action classes

 In the RE04ActionInitialization class, three user action classes are commented out,
 i.e.
- RE04EventAction
- RE04TrackingAction
- RE04SteppingAction

By using RE04SteppingAction, you will
see a material name which a track sees for each step.
By using RE04EventAction and RE04TrackingAction, you will
see the similar information for all trajectories of one
event.

