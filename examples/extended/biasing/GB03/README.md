\page ExampleGB03 Example GB03

## Geometry based biasing

This example illustrates a use of generic biasing classes to implement a
technique near to "geometry importance biasing".

The geometry is the same than in EM tests, with the sampling calorimeter
made of a series of layers of absorber and gap.

The biasing applies to neutrons only.

Instead of explicitly assigning "importance" values to the layers, we
split neutrons moving forward and kill the ones moving backward, when they
reach the exit of an absorber volume.

The splitting factor can be controlled by command line, eg:
```
/GB03/biasing/setSplittingFactor 2
```

which also determines the killing probability : 1/(splitting factor).

It can be seen than when defining 10 layers (see exampleGB03.in), a
splitting factor 2 works fine : we don't suffer from under- or over-splitting.
If going to 20 layers, then a splitting with a factor 2 is too large,
and the biasing suffers from over-splitting. (And we can not go lower than
"2", which would mean "1" and hence, no biasing...)

To alleviate the over-splitting, we introduce a probability to apply the
splitting (and killing) (this is one solution, others can be considered), that
can be changed as:

```
/GB03/biasing/setApplyProbability 0.5
```

With above value, we can see that we recover a satisfactory biasing scheme,
with neutrons penetrating the entire setup, without over-splitting.

The commands in ```/GB03/biasing``` section are available only if application
started with biasing "on".

The classes involved are:

- GB03BOptnSplitOrKillOnBoundary : which is the biasing operation making
  the splitting and killing;
- GB03BOptrGeometryBasedBiasing : which is the biasing operator, making
  decision to use above operation, and configuring it, passing it the
  splitting factor and probability to apply the biasing.

## HOW TO START ?

The example can be executed in "batch" mode if macro file name is specified
or interactive mode.

To run it:
```
  ./exampleGB03 [-m macro_file] [-b [on|off]]
```

Without parameters the application runs in interactive mode with biasing "on".

## THE OUTPUT

During the run the histograms of energy distributions for neutrons and gammas
after the "shield" volume and their positions along x axis are filed and
saved to GB03.root file.

Depending on ```/GB03/verbose``` flag

For *verbose >= 1* all particles penetrating the shield were counted and
printed with their weights and Emin Emax values.

The track #, particle, kinetic energy and position on each step are printed
if *verbose >= 2*.

The particle count per thread is also printed if *verbose >= 3*.
