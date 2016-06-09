Geant4 MPI Interface
====================

Author:
Koichi Murakami (KEK) / Koichi.Murakami@kek.jp


About the interface
===================
G4MPI is a native interface with MPI libraries. The directory contains 
a Geant4 UI library and a couple of parallelized examples.
Using this interface, users applications can be parallelized with
different MPI compliant libraries, such as OpenMPI, LAM/MPI, MPICH2 and so on.

System Requirements:
--------------------

### MPI Library

The MPI interface can work with MPI-compliant libraries, 
such as Open MPI, LAM/MPI, Intel MPI etc.

For example, the information about Open MPI can be obtained from
http://www.open-mpi.org/

### CMake

CMake is used to build G4MPI library, that co-works with Geant4 build system.

### Optional

ROOT for histogramming/analysis

- - -

How to build G4MPI
==================
To build G4MPI library, use CMake on Geant4 library installed with CMake build.

Check `CMakeList.txt`, especially the following two variables 
should be taken care to match your MPI library

    eg.  
    set(CMAKE_CXX_COMPILER mpicxx)  
    #set(CMAKE_CXX_INCLUDE_PATH ) set if necessary

Follow these commands,

    > cd source
    > mkdir build
    > cd build
    > cmake -DGeant4_DIR=<your Geant4 install path>/lib64/Geant4-V.m.n ..
      (V.m.n is the version of Geant4, eg. Geant4-9.6.0)
    > make
    > make install

The library and header files will be installed on the installation directory
of Geant4.

- - - 

How to use
==========

How to make parallel applications
---------------------------------

An example of a main program:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cc}
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
   `lamboot` (LAM) / `mpdboot` (MPICH2).

How to run
----------
For example,

    > mpiexec -n # <your application>

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

- - -

Examples
========
There are a couple of examples for Geant4 MPI applications.

In some cases, you need to set some additional environment variables
for running examples:

- *G4LEDATA* : directory path for low energy EM data
- *G4LEVELGAMMADATA* : directory path for photon evapolation
- *G4SAIDXSDATA* : directory path for nucleon cross section data


For using ROOT libraries

- *ROOTSYS* : root path of the ROOT package

exMPI01
-------
A simple application.

**Configuration:**

- Geometry : chamber / calorimeter
- Primary : particle gun (200 MeV electron as default)
- Physics List : standard EM

**Features:**
- Particles are transported in a geometry without any scoring.
- Learn how to parallelized your G4 session.

exMPI02 (ROOT application)
--------------------------
An example of dosimetry in a water phantom.

**Configuration:**
- Geometry     : water phantom
- Primary      : broad beam (200 MeV proton)
- Physics List : FTFP_BERT
- Analysis     : ROOT histogramming

**Features:**
- Score dose distribution in a water phantom.
- Learn how to parallelized your applications.
- Create a ROOT file containing histograms/trees in each node.
