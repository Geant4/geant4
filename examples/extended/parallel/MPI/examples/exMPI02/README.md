MPI/Examples : exMPI02
======================

Description
-----------
An example of dosimetry in a water phantom.

### Configuration:
- Geometry     : water phantom
- Primary      : broad beam (200 MeV proton)
- Physics List : FTFP_BERT
- Analysis     : ROOT histogramming

The environment variables, *G4LEDATA*, *G4LEVELGAMMADATA* and *G4SAIDXSDATA*
are required for data files.

The enviromnet variable, *ROOTSYS*, is set to the root path of the ROOT package.

### Features:
- Score dose distribution in a water phantom.
- Learn how to paralleized your applications.
- Create a ROOT file containing histograms/trees in each node.
  Each slave node generate a ROOT file, whose file name is different 
  from each other.


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
