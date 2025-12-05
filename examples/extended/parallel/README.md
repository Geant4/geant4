\page Examples_parallel Category "parallel"

This directory includes example applications to demonstrate the usage of
different techniques for achieving event parallelism with Geant4.

- \ref Examples_MPI

  The directory contains a native interface with MPI libraries, 
  as a Geant4 UI library, and a set of parallelized examples.
  Using this interface, users applications can be parllelized with 
  different MPI compliant libraries, such as LAM/MPI, MPICH2, OpenMPI,
  and so on. 

- \ref ExampleThreadsafeScorers

 This example demonstrates a very simple application where an energy
 deposit and # of steps is accounted in thread-local (i.e. one instance per
 thread) hits maps with underlying types of plain-old data (POD) and global
 (i.e. one instance) hits maps with underlying types of atomics.
