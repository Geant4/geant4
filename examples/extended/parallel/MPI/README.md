Geant4 MPI Interface
====================

Author:
Koichi Murakami (KEK) / Koichi.Murakami@kek.jp
Andrea Dotti (SLAC) / adotti@slac.stanford.edu


About the interface
===================
G4MPI is a native interface with MPI libraries. The directory contains 
a Geant4 UI library and a couple of parallelized examples.
Using this interface, users applications can be parallelized with
different MPI compliant libraries, such as OpenMPI, MPICH2 and so on.

System Requirements:
--------------------

### MPI Library

The MPI interface can work with MPI-compliant libraries, 
such as Open MPI, MPICH, Intel MPI etc.

For example, the information about Open MPI can be obtained from
http://www.open-mpi.org/

MPI support:
------------
G4mpi has been tested with the following MPI flavors:
 *   OpenMPI 1.8.1
 *   MPICH 3.2
 *   Intel MPI 5.0.1
 
### CMake

CMake is used to build G4MPI library, that co-works with Geant4 build system.

### Optional (for exMPI02)

ROOT for histogramming/analysis

How to build G4MPI
==================
To build G4MPI library, use CMake on Geant4 library installed with CMake build.

Follow these commands,

    > mkdir build
    > cd build
    > cmake -DGeant4_DIR=<your Geant4 install path>/lib[64]/Geant4-V.m.n \
     -DCMAKE_INSTALL_PREFIX=<where-G4mpi-lib-will-be-installed> \
     <g4source>/examples/extended/parallel/MPI/source
    > make
    > make install

The cmake step will try to guess where MPI is installed, mpi executables should 
be in PATH. You can specify CXX and CC environment variables to your specific 
mpi wrappers if needed.

The library and header files will be installed on the installation directory 
specified in CMAKE_INSTALL_PREFIX
a CMake configuration file will also be installed 
(see examples on how to compile an application using G4mpi)

How to use
==========

How to make parallel applications
---------------------------------

An example of a main program:

```c++
#include "G4MPImanager.hh"
#include "G4MPIsession.hh"

int main(int argc,char** argv)
{
  // At first, G4MPImanager/G4MPIsession should be created.
  G4MPImanager* g4MPI= new G4MPImanager(argc,argv);
      
  // MPI session (G4MPIsession) instead of G4UIterminal
  G4MPIsession* session= g4MPI-> GetMPIsession();
      
  // user application setting
  G4RunManager* runManager= new G4RunManager();

  ....

  // After user application setting, just start a MPI session.
  MPIsession treats both interactive and batch modes.
  session-> SessionStart();

  // Finally, terminate the program
  delete g4MPI;
  delete runManager;
}
```

How to compile
---------------

Using cmake, assuming G4mpi library is installed in path _g4mpi-path_ and 
Geant4 is installed in _g4-path_:

    > mkdir build
    > cd build
    > cmake -DGeant4_DIR=<g4-path>/lib[64]/Geant4-V.m.n \
      -DG4mpi_DIR=<g4mpi-path>/lib[64]/G4mpi-V.m.n \
      <source>
    > make

Check provided examples: under examples/extended/parallel/MPI/examples for an
example of CMakeLists.txt file to be used.

### Notes about session shell

LAM/MPI users can use "G4tcsh" as an interactive session shell.
For other users (Open MPI/MPICH2), plesae use G4csh (default).

In case of OpenMPI, *LD_LIBRARY_PATH* for OpenMPI runtime libraries
should be set at run time. Alternatively, you can add this path
to the dynamic linker configuration using `ldconfig`.
(needs sys-admin authorization)


MPI runtime Environment
-----------------------
1. Make hosts/cluster configuration of your MPI environment.
2. Launch MPI runtime environment, typically executing
   `lamboot` (LAM) / `mpdboot` (MPICH2) / `mpd` (Intel).

How to run
----------
For example,

    > mpiexec -n # <your application>
<p>Replace mpicxx with your MPI compiler wrapper if you need to specify which one to use.</p>

Instead, `mpirun` command is more convenient for LAM users.


MPI G4UI commands
-----------------
G4UI commands handling the G4MPI interface are placed in /mpi/.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Command directory path : /mpi/
    
Guidance :
MPI control commands
    
 Sub-directories :
 Commands :
   verbose *        Set verbose level.
   status *         Show mpi status.
   execute *        Execute a macro file. (=/control/execute)
   beamOn *         Start a parallel run w/ thread.
   .beamOn *        Start a parallel run w/o thread.
   masterWeight *   Set weight for master node.
   showSeeds *      Show seeds of MPI nodes.
   setMasterSeed *  Set a master seed for the seed generator.
   setSeed *        Set a seed for a specified node.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Notes:
While "/run/beamOn" and "/mpi/beamOn" commands invoke beam-on in background,
so you can input UI commands even while event processing. Note that drawing
tracks in OpenGL with these commands causes a crash. Please use /mpi/.beamOn
command instead.

The original "/control/execute" and "/run/beamOn" are overwritten
with "/mpi/execute" and "/mpi/beamOn" commands respectively,
that are customized for the MPI interface.

Examples
========
There are a couple of examples for Geant4 MPI applications.

For using ROOT libraries (exMPI02)

- *ROOTSYS* : root path of the ROOT package

exMPI01
-------
A simple application.

**Configuration:**

- Geometry : chamber / calorimeter
- Primary : particle gun (200 MeV electron as default)
- Physics List : FTFP_BERT

**Features:**
- Particles are transported in a geometry without any scoring.
- Learn how to parallelized your G4 session.

exMPI02 (ROOT application)
--------------------------
An example of dosimetry in a water phantom.
Note: due to limited MT support in ROOT, in this example
      MT is disabled, but the code is migrated to MT, ready
      for MT when ROOT will support MT.
      For an example of MT+MPI take a look at exMPI03

**Configuration:**
- Geometry     : water phantom
- Primary      : broad beam (200 MeV proton)
- Physics List : FTFP_BERT
- Analysis     : ROOT histogramming

**Features:**
- Score dose distribution in a water phantom.
- Learn how to parallelized your applications.
- Create a ROOT file containing histograms/trees in each node.

exMPI03 (merging of histograms via MPI)
---------------------------------------
This example is the same as exMPI02 with the following
differences:
- It uses Geant4 analysis instead of ROOT for histogramming
- It shows how to merge, using g4tools, histograms via MPI
  so that the entire statistics is accumulated in a single output file
- It also shows how to merge G4Run objects from different ranks and 
  how to merge scorers
- MT is enabled.
- Root output files from application run with `mpiexec -n 3`
  - dose-merged.root - merged histograms
  - dose-rank0,1,2 - histograms data collected on rank 0, 1,2 before merge

exMPI04 (merging of ntuples via MPI)
---------------------------------------
This example is the same as exMPI03 with added ntuple.
- It uses Geant4 analysis for histogramming and ntuples.
- It shows how to merge, using g4tools, ntuples via MPI in sequential mode,
  so that the entire statistics is accumulated in a single output file.
- If MT is enabled, the ntuples are merged from threads to 
  files per ranks.
- Combined MT + MPI merging is not yet supported.
- Merging ntuples is actually supported only with Root output format.

- Root output files from application run with `mpiexec -n 4`
  - Sequential application:
    (3 working ranks, the last rank dedicated for collecting ntuple data)
    - dose-merged.root - merged histograms from ranks 0, 1 and 2
    - dose-rank1,2.root - histograms data collected on rank N before merge
    - dose-rank3  - ntuples merged from ranks 0, 1 and 2
  - MT application:
    (4 working ranks)
    - dose-merged.root - merged histograms from ranks 0, 1, 2 and 3
    - dose-rank0,1,2,3.root -  histograms data collected on rank N before merge;
         ntuples merged on each rank N from rank threads

