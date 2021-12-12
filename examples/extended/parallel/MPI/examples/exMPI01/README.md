MPI/Examples : exMPI01
======================

Description
-----------
A simple application

### Configuration:

- Geometry     : chamber / calorimeter
- Primary      : particle gun (200 MeV electron as default)
- Physics List : FTFP_BERT

### Features:

- Particles are transported in a geometry without any scoring.
- Learn how to parallelized your G4 session.


How to build
------------
Use CMake on Geant4 library installed with CMake build.

This example requires G4mpi library to be installed
(see examples/extended/parallel/MPI/source/REDME.md)

Follow these commands,

    > mkdir build
    > cd build
    > cmake -DG4mpi_DIR=<where-G4mpi-wasintalled>/lib[64]/G4mpi -DCMAKE_CXX_COMPILER=mpicxx \
      -DGeant4_DIR=<your Geant4 install path>/lib[64]/Geant4-V.m.n <path-to-source>
      (V.m.n is the version of Geant4, eg. Geant4-9.6.0)
    > make
    > make install

Replace mpicxx with your MPI-compiler wrapper.
