
     ================================================================================================
                    Geant4 - an Object-Oriented Toolkit for Simulation in HEP
     ================================================================================================

                                        dsbandrepair
                                         ---------

                        **A Geant4-DNA application for simulating early DNA damage**

# AUTHORS
L. T. Anh, Y. Perrot, C. Villagrasa, S. Meylan, H. N. Tran

(\*) contact: yann.perrot@irsn.fr or carmen.villagrasa@irsn.fr

# Introduction

“dsbandrepair” is a Geant4-DNA simulation chain for evaluating the early radiation-induced DNA damage.
The first development of the simulation chain was carried out by Meylan et al. in 2017 (Sci. Rep. 2017 7:11923)
The "extended/medical/dna/dnadamage1" example is a simplified version of "dsbanrepair"

“dsbandrepair” supports all types of DNA geometries constructed with DNAFabric (Comput. Phys. Comm. 2016 204:159-169).
Geometries for human cell nuclei (fibroblast, endothelium) and yeast were provided along with the release of “dsbandrepair”.
Users can use a free version of DNAFabric (https://bitbucket.org/sylMeylan/opendnafabric/src/master/) to create customed geometries. Or they can contact Y. Perrot for specific geometries.
The geometric models are constructed from 10 voxels to form a continuous chromatin fiber for each chromosme including heterochromatin (VoxelStraight, VoxelRight,...) and euchromatin (VoxelStraight2, VoxelRight2,...) distribution (Med. Phys. 2019 46:1501-1511).

Physical stage and chemical stage allow the calculation of direct and indirect Strand Breaks in the whole nucleus.

Furthermore, the Two Lesion Kinetic model (Radiat. Res. 2001 156:365-378) and the Local Effect Model IV (Radiat. Res. 2013 180:524-538) were also included to allow users calculate the survival fraction and un-rejoined DSBs.
The Belov's model (J. Theo. Biol. 2015 366:115-130) for double-strand breaks repair is provided but has not been compared to experimental data.
 

# How to build and run

To build dsbandrepair, in the terminal, use:
*  shell$ mkdir build
*  shell$ cd build
*  shell$ cmake /path-to/dsbandrepair  
(Or if users don't want to download geometry files while compiling the dsbandrepair, use: cmake -DDOWNLOAD_GEOMETRY=FALSE /path-to/dsbandrepair )
*  shell$ make  (or 'make -jN' with N = 1,2,3 .... )

And to run:
*  shell$ ./dsbandrepair dsbandrepair.in

where dsbandrepair.in is a macrofile. User can change it to his/her own macrofile.

Note that: dsbandrepair was designed in a modular way that offers users to run physical stage chemcal stage independently. By default, dsbandrepair runs in physical stage mode. To run chemical stage, use :
*  shell$ rm -rf chem_ouput
*  shell$ ./dsbandrepair chem.in chem 

where chem.in is a macrofile. User can change it to his/her own macrofile.

## Running with mpi library

To improve the simulation in term of computational time, user can run dsbandrepair with mpi library.

MPI interface: Thanks to the work of K. Murakami and A. Dotti (DOI: https://doi.org/10.1109/NSSMIC.2015.7581867), an MPI interface was introduced into Geant4 and it’s now used in this work (see "/examples/extended/parallel/MPI"). User has to follow this example to install g4mpi library.

To compile the "dsbandrepair" with g4mpi:
*  shell$ mkdir build
*  shell$ cd build
*  shell$ cmake -DUSE_MPI=TRUE -DG4mpi_DIR=<g4mpi-path>/lib[64]/G4mpi-V.m.n /path-to/dsbandrepair  
*  shell$ make  (or 'make -jN' with N = 1,2,3 .... )
And to run:
*  shell$ mpiexec -np $nranks ./dsbandrepair dsbandrepair.in

where $nranks is the number of mpi processes you want to run.

Or ro run chemical stage:

*  shell$ rm -rf chem_ouput
*  shell$ mpiexec -np $nranks ./dsbandrepair chem.in chem 

# Analyzing results
To run "analysis" module, in the "build" directory, build this module with the commands:
* shell$  mkdir analysis
* cd analysis
* cmake /path/to/analysis
* make 
* cd ../

At this point, user can launch the analysis module:
* shell$  ./analysis/runAna  
or
* shell$  ./analysis/runAna macrofile

where the macro file allows user to interact with the code. 
Example: ./analysis/runAna analysis.in 

## Outputs
By default, the output of "Analysis" module will be written in 4 different text files:
* SB results: this text file contains all SB results, such as total SB, direct and indirect SBs, SSB and DSB.
* SDD format: All damages are written in SDD format (Radiat. Res. 2019 191:11). File name starts with "SDD_" 
* TLK result: File name starts with "TLK_". This file contains results from TLK model.
* LEM-IV result: File name starts with "LEMIV_". This file contains results from LEMIV model.

# Maro files:
Some macro files are provided along with this code, user can change them based on their own needs.

* macro files for physical stage:
    * dsbansrepair.in : This macro is for a light geometry for testing the code
    * fibroblast.in: This macro is for fibroblast cell nucleus.
    * endophys.in: This macro is for endothelium cell nucleus.
    * yeastphys.in: This macro is for yeast cell nucleus.
* macro files for chem stage:
    * chem.in
* macro files for analysis module:
    * analysis.in: allows user to set parameter for scoring, classifying damages and setting repair models parameters.

# An alternative example for DNA damage calculation can be found in examples/extended/medical/dna/moleculardna

# Acknowledgments
The transition from the initial simulation chain of Meylan et al. to a version adapted for a Geant4 example benefited from funds from the BioRad3 project financed by the ESA (grant DAR 4000132935/21/NL/CRS)
